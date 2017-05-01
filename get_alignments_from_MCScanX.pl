#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;
use Bio::AlignIO;
use Sort::Naturally;
use File::Path 'rmtree';
use Bio::Align::Utilities qw(aa_to_dna_aln);

my $usage = "
SYNOPSIS
  Generates codon alignments for pairs of genes with a given minimum Ks from an MCScanX collinearity file (annotated with Ka/Ks).

OPTIONS
  -i|--in      [FILE] : MCScanX collinearity file
  -p|--prot    [FILE] : fasta file of protein sequences
  -c|--cds     [FILE] : fasta file of corresponding CDS (nucleotide)
  -o|--out      [STR] : outfile prefix (default = 'ALN'); alignments will be written to ALN.<NUM>.fasta
  -d|--outdir   [DIR] : dirname to save alignments (default = 'alignments')
  -k|--minks  [FLOAT] : minimum Ks between genes (default >= 0.5)
  -t| --threads [INT] : number of clutalo threads (default = 1)
  -h|--help           : this message

USAGE

\n";

my ($infile, $proteinfile, $cdsfile, $help);
my $minKs = 0.5;
my $out = "ALN";
my $outdir = "alignments";
my $threads = 1;

GetOptions (
  'in|i=s'      => \$infile,
  'prot|p=s'    => \$proteinfile,
  'cds|c=s'     => \$cdsfile,
  'o|out:s'     => \$out,
  'd|outdir'    => \$outdir,
  'minks|k:f'   => \$minKs,
  'threads|t:i' => \$threads,
  'help|h'      => \$help
);

die $usage unless ($infile && $proteinfile && $cdsfile);
die $usage if $help;
print STDERR "[INFO] Collinearity file: $infile\n";
print STDERR "[INFO] Minimum Ks: $minKs\n";
print STDERR "[INFO] Number of clustalo threads set to: $threads\n";
print STDERR "[INFO] Parsing collinearity file...\n";

## parse collinearity file:
open (my $IN, $infile) or die "Cannot open $infile: $!\n\n";
my %pairs;
while (<$IN>) {
  chomp;
  next if /^\#/;
  $_ =~ s/^\s+|\s+$//g; ##remove leading and trailing whitespaces
  my @F = split (/\s+/,$_);
  ## work from end of array as columns in collinearity file not consistently formatted
  $pairs{$F[-5]}=$F[-4] if $F[-1] >= $minKs;
}
close $IN;
print STDERR "[INFO] Number of pairs with Ks >= $minKs: ".commify(scalar(keys(%pairs)))."\n";

## parse proteins and CDS:
my (%protein_hash, %cds_hash);
my $in_p = Bio::SeqIO->new( -file => $proteinfile, -format => 'fasta' );
while (my $seq = $in_p->next_seq() ) {
  $protein_hash{$seq->display_id()} = $seq->seq();
}
print "[INFO] Fetched ".commify(scalar(keys %protein_hash))." proteins from $proteinfile\n";
my $in_c = Bio::SeqIO->new( -file => $cdsfile, -format => 'fasta' );
while (my $seq = $in_c->next_seq() ) {
  $cds_hash{$seq->display_id()} = $seq->seq();
}
print "[INFO] Fetched ".commify(scalar(keys %cds_hash))." CDS from $cdsfile\n";
die "[ERROR] No sequences found in $proteinfile or $cdsfile!\n" if ((scalar(keys %protein_hash) == 0) || (scalar(keys %cds_hash) == 0));

## generate alignments for all pairs in %pairs:
my $n = 1;
open (my $OUT, ">$infile.map") or die "[ERROR] Cannot open $infile.map: $!\n\n";
print $OUT "NUM\tSEQ1\tSEQ2\n";
foreach (nsort keys %pairs) {
  print STDERR "\r[INFO] Working on pair \#$n: $_, $pairs{$_}";$| = 1;
  print $OUT "$n\t$_\t$pairs{$_}\n";

  ## fetch proteins and print to temp file
  open (my $PRO, ">clustal.pro") or die $!;
  if ((exists($protein_hash{$_})) && (exists($protein_hash{$pairs{$_}}))) {
    print $PRO ">$_\n$protein_hash{$_}\n>$pairs{$_}\n$protein_hash{$pairs{$_}}";
    close $PRO;
  } else {
    die "[ERROR] Protein ID '$_' or '$pairs{$_}' not found in file $proteinfile!\n";
  }

  ## make CDS hash of nucleotides
  my %cds_seqs;
  if ((exists($cds_hash{$_})) && (exists($cds_hash{$pairs{$_}}))) {
    $cds_seqs{"$_"} = Bio::Seq->new( -display_id => "$_", -seq => $cds_hash{$_} );
    $cds_seqs{"$pairs{$_}"} = Bio::Seq->new( -display_id => "$pairs{$_}", -seq => $cds_hash{$pairs{$_}} );
  } else {
    die "[ERROR] CDS ID '$_' or '$pairs{$_}' not found in file $cdsfile!\n";
  }

  ## run alignment
  if (system ("clustalo --infile=clustal.pro --outfile=clustal.aln --force --threads=$threads") != 0) { die "[ERROR] Problem with clustalo!\n"; }

  ## fetch alignment, backtranslate to nucleotides & write
  my $get_prot_aln = Bio::AlignIO -> new(-file=>"clustal.aln", -format=>"fasta");
  my $prot_aln = $get_prot_aln -> next_aln();
  my $dna_aln = aa_to_dna_aln($prot_aln, \%cds_seqs);
  my $write_dna_aln = Bio::AlignIO -> new(-file=>">$out.$n.fasta", -format=>"fasta");
  $write_dna_aln->write_aln($dna_aln);
  $n++;
}
close $OUT;

## file cleanup:
if (-d $outdir) {
  rmtree([ "$outdir" ]);
  mkdir $outdir;
  system("mv $out* $outdir");
} else {
  mkdir $outdir;
  system("mv $out* $outdir");
}
system ("rm clustal*");
print "\n[INFO] Finished on ".`date`."\n";

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
