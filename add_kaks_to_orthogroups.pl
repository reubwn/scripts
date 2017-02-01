#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq::EncodedSeq;
use Bio::Align::DNAStatistics;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Algorithm::Combinatorics qw(combinations);

my $usage = "
SYNOPSIS:
  Calculates Ka and Ks for all pairs of genes from an OrthoFinder Orthogroups.txt file.

OPTIONS:
  -i|--orthogroups [FILE]   : Orthogroups.txt file from OrthoFinder
  -1|--sp1         [STRING] : species 1 identifier to select genes from (genome ID in protein identifier)
  -2|--sp2         [STRING] : species 2 identifier to select genes from (optional) ** NOT IMPLEMENTED YET! **
  -p|--protein     [FILE]   : fasta file of protein sequences
  -c|--cds         [FILE]   : fasta file of corresponding CDSs (nucleotide)
  -o|--out         [FILE]   : output file to write to (default: 'inputfilename.kaks')
  -n|--noheader             : don't print header to outfile (default: do print it)
  -h|--help                 : print this message

USAGE:
  Generate Ka and Ks for all pairs of genes beginning with ID 'ARIC' from file 'Orthogroups.txt':
  >> add_ka_and_ks_to_orthogroups.pl -i Orthogroups.txt -1 ARIC -p <(cat proteins/*) -c <(cat cds/*) -o ARIC.kaks
\n";

my ($orthogroups, $sp1, $sp2, $proteinfile, $cdsfile, $outfile, $noheader, $help);

GetOptions (
  'orthogroups|i=s' => \$orthogroups,
  'sp1|1=s'         => \$sp1,
  'sp2|2:s'         => \$sp2,
  'protein|p=s'     => \$proteinfile,
  'cds|c=s'         => \$cdsfile,
  'outfile|o:s'     => \$outfile,
  'noheader|n'      => \$noheader,
  'help|h'          => \$help
);

die $usage if $help;
die $usage unless ($orthogroups && $sp1 && $proteinfile && $cdsfile);

$outfile = $orthogroups."kaks" unless ($outfile);

## parse proteins and CDSs
my (%protein_hash, %cds_hash);
my $in_p = Bio::SeqIO->new( -file => $proteinfile, -format => 'fasta' );
while (my $seq = $in_p->next_seq() ) {
  $protein_hash{$seq->display_id()} = $seq->seq();
}
print "[INFO] Fetched ".(scalar(keys %protein_hash))." protein seqs from $proteinfile\n";
my $in_c = Bio::SeqIO->new( -file => $cdsfile, -format => 'fasta' );
while (my $seq = $in_c->next_seq() ) {
  $cds_hash{$seq->display_id()} = $seq->seq();
}
print "[INFO] Fetched ".(scalar(keys %cds_hash))." CDS seqs from $cdsfile\n";
die "[ERROR] No sequences found in $proteinfile or $cdsfile!\n" if ((scalar(keys %protein_hash) == 0) || (scalar(keys %cds_hash) == 0));

open (my $OUT, ">$outfile") or die $!;
open (my $GROUPS, $orthogroups) or die $!;
print $OUT "OG\tgene1\tgene2\tKa\tKs\n" unless ($noheader);
while (<$GROUPS>) {
  my @F = split (m/\s+/);
  my $OG_num = shift @F;
  my @seqs_of_interest;
  foreach (@F) {
    if ($_ =~ m/^$sp1/) { ## push if $sp1 matches START of gene ID
      push @seqs_of_interest, $_;
    }
  }

  ## get all pairs
  my $iter = combinations(\@seqs_of_interest, 2);

  ## iterate thru pairs (NB $combo is an arrayref)
  while (my $combo = $iter->next()) {
    print "[INFO] Working on pair: @$combo\n";

    ## fetch proteins and print to temp file
    open (my $PRO, ">clustal.pro") or die $!;
    if (($protein_hash{@{$combo}[0]}) && ($protein_hash{@{$combo}[1]})) {
      print $PRO ">@{$combo}[0]\n$protein_hash{@{$combo}[0]}\n>@{$combo}[1]\n$protein_hash{@{$combo}[1]}";
      close $PRO;
    } else {
      die "[ERROR] Protein ID '@{$combo}[0]' or '@{$combo}[1]' not found in file $proteinfile!\n";
    }

    ## make CDS hash of nucleotides
    my %cds_seqs;
    if (($cds_hash{@{$combo}[0]}) && ($cds_hash{@{$combo}[1]})) {
      $cds_seqs{"@{$combo}[0]"} = Bio::Seq->new( -display_id => "@{$combo}[0]", -seq => $cds_hash{@{$combo}[0]} );
      $cds_seqs{"@{$combo}[1]"} = Bio::Seq->new( -display_id => "@{$combo}[1]", -seq => $cds_hash{@{$combo}[1]} );
    } else {
      die "[ERROR] CDS ID '@{$combo}[0]' or '@{$combo}[1]' not found in file $cdsfile!\n";
    }

    ## run alignment
    if (system ("clustalo --infile=clustal.pro --outfile=clustal.aln --force") != 0) { die "[ERROR] Problem with clustalo!\n"; }

    ## fetch alignment and backtranslate to nucleotides
    my $get_prot_aln = Bio::AlignIO -> new(-file=>'clustal.aln', -format=>'fasta');
    my $prot_aln = $get_prot_aln -> next_aln();
    my $dna_aln = aa_to_dna_aln($prot_aln, \%cds_seqs);

    ## get Ka (Dn), Ks (Ds) values
    eval {
      my $stats = Bio::Align::DNAStatistics->new();
      my $result = $stats->calc_all_KaKs_pairs($dna_aln);
      my ($Da, $Ds, $Dn, $N, $S, $S_d, $N_d);
      for my $an (@$result) {
        for (sort keys %$an ) {
          next if /Seq/;
          if($_ eq "D_n"){$Dn = $an->{$_}};
          if($_ eq "D_s"){$Ds = $an->{$_}};
          if($_ eq "S_d"){$S_d = $an->{$_};}
          if($_ eq "N_d"){$N_d = $an->{$_};}
          if($_ eq "S"){$S = $an->{$_};}
          if($_ eq "N"){$N = $an->{$_};}
        }
      }
      $Dn = -2 unless ($Dn); ## default values
      $Ds = -2 unless ($Ds);
      #if($Dn !~ /\d/){$Dn = -2;}
      #if($Ds !~ /\d/){$Ds = -2;}

      ## print to file
      print $OUT "$OG_num\t@{$combo}[0]\t@{$combo}[1]\t$Dn\t$Ds\n";
    }
  }
}
close $GROUPS;
close $OUT;
system ("rm clustal*");

__END__
