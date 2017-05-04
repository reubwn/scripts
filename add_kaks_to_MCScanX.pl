#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq::EncodedSeq;
use Bio::Align::DNAStatistics;
use Bio::Align::Utilities qw(aa_to_dna_aln);

local $SIG{__WARN__} = sub { warn $_[0] unless $_[0] =~ /WARNING/}; ##suppress WARNINGs from BioPerl

my $usage = "
SYNOPSIS:
  Calculates Ka and Ks for pairs of genes from an MCScanX collinearity file.

OPTIONS
  -i|--in      [FILE] : MCScanX collinearity file
  -p|--prot    [FILE] : fasta file of protein sequences
  -c|--cds     [FILE] : fasta file of corresponding CDS (nucleotide)
  -t|--threads  [INT] : number of aligner threads (default = 1)
  -h|--help           : this message

USAGE
\n";

my ($infile, $proteinfile, $cdsfile, $help);
my $threads = 1;

GetOptions (
  'in|i=s'      => \$infile,
  'prot|p:s'    => \$proteinfile,
  'cds|c:s'     => \$cdsfile,
#  'o|out:s'     => \$outfile,
  'threads|t:i' => \$threads,
  'help|h'      => \$help
);

die $usage if $help;
die $usage unless ($infile && $proteinfile && $cdsfile );
my $outfile = "$infile.kaks";

print STDERR "[INFO] Collinearity file: $infile\n";
print STDERR "[INFO] Proteins file: $proteinfile\n";
print STDERR "[INFO] CDS file: $cdsfile\n";
print STDERR "[INFO] Parsing protein and CDS files...\n";

## parse proteins and CDS:
my (%protein_hash, %cds_hash);
my $in_p = Bio::SeqIO->new( -file => $proteinfile, -format => 'fasta' );
while (my $seq = $in_p->next_seq() ) {
  $protein_hash{$seq->display_id()} = $seq->seq();
}
die "[ERROR] No sequences found in $proteinfile!\n" if (scalar(keys %protein_hash) == 0);
print "[INFO] Fetched ".commify(scalar(keys %protein_hash))." proteins from $proteinfile\n";
my $in_c = Bio::SeqIO->new( -file => $cdsfile, -format => 'fasta' );
while (my $seq = $in_c->next_seq() ) {
  $cds_hash{$seq->display_id()} = $seq->seq();
}
die "[ERROR] No sequences found in $cdsfile!\n" if (scalar(keys %cds_hash) == 0);
print "[INFO] Fetched ".commify(scalar(keys %cds_hash))." CDS from $cdsfile\n";

print STDERR "[INFO] Parsing collinearity file...\n";

## parse collinearity file:
open (my $OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
open (my $IN, $infile) or die "Cannot open $infile: $!\n\n";
my ($trimmed_seqs,$na) = (0,0);
my $n = 1;
while (<$IN>) {
  if ($_ =~ m/^\#/) {
    print $OUT $_;
  } else {
    chomp;
    my $LINE = $_;
    $LINE =~ s/^\s+|\s+$//g; ##remove leading and trailing whitespaces
    my @F = split (/\s+/, $LINE);
    #print STDERR "$LINE\t$F[-3]\t$F[-2]\n";
    my $gene1 = $F[-3]; ##work from end of array as columns in collinearity file not consistently formatted
    my $gene2 = $F[-2];

    print STDERR "\r[INFO] Working on pair \#$n: $gene1, $gene2";$| = 1;

    ## fetch proteins and print to temp file:
    open (my $PRO, ">temp.faa") or die $!;
    if ((exists($protein_hash{$gene1})) && (exists($protein_hash{$gene2}))) {
      print $PRO ">$gene1\n$protein_hash{$gene1}\n>$gene2\n$protein_hash{$gene2}";
      close $PRO;
    } else {
      die "[ERROR] Protein ID '$gene1' or '$gene2' not found in file $proteinfile!\n";
    }

    ## run clustalo alignment:
    if (system ("clustalo --infile=temp.faa --outfile=temp.aln --force --threads=$threads") != 0) { die "[ERROR] Problem with clustalo!\n"; }

    ## make CDS hash of nucleotides:
    my %cds_seqs;
    if ((exists($cds_hash{$gene1})) && (exists($cds_hash{$gene2}))) {
      $cds_seqs{"$gene1"} = Bio::Seq->new( -display_id => "$gene1", -seq => $cds_hash{$gene1} );
      $cds_seqs{"$gene2"} = Bio::Seq->new( -display_id => "$gene2", -seq => $cds_hash{$gene2} );
    } else {
      die "[ERROR] CDS '$gene1' or '$gene2' not found in file $cdsfile!\n";
    }
    ## check if all CDS seqs are multiple of 3:
    foreach (keys %cds_seqs) {
      if ($cds_seqs{$_}->length() % 3 != 0) {
        print STDERR "\n[WARN] Seq ".$cds_seqs{$_}->display_id()." is not a multiple of 3; will trim ".($cds_seqs{$_}->length() % 3)." bases from 3' end\n";
        my $trimmed = $cds_seqs{$_}->subseq(1,(($cds_seqs{$_}->length()) - ($cds_seqs{$_}->length() % 3))); ##trims remainder off 3' end; returns a STRING, annoyingly
        $cds_seqs{$_} = Bio::Seq->new( -display_id => $_, -seq => $trimmed ); ##replace old seq with trimmed seq
        $trimmed_seqs++;
      }
    }

    ## fetch alignment, backtranslate to nucleotides:
    my $get_prot_aln = Bio::AlignIO -> new(-file=>"temp.aln", -format=>"fasta");
    my $prot_aln = $get_prot_aln -> next_aln();
    my $dna_aln = aa_to_dna_aln($prot_aln, \%cds_seqs);

    ## get Ka (Dn), Ks (Ds) values:
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
      $Dn = -2 unless ($Dn); ##default values
      $Ds = -2 unless ($Ds);
      $na++ if ( ($Dn == -2) || ($Ds == -2) );

      ## annotate input file
      print $OUT "$_\t$Dn\t$Ds\n";
    };
    $n++;
  }
}
close $IN;
close $OUT;
system ("rm temp.*"); ##remove last temp files.

print STDERR "\n";
print STDERR "[INFO] Total number of pairs: $n\n";
print STDERR "[INFO] Number of seqs trimmed because % 3 != 0: $trimmed_seqs\n";
print STDERR "[INFO] Number of pairs for which Ka or Ks was not calculated: $na\n";
print STDERR "[INFO] Finished on ".`date`."\n";

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
