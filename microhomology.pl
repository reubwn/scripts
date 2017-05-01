#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Sort::Naturally;

use Bio::Seq;
use Bio::AlignIO;

my $usage = "
SYNOPSIS
  Calculates the number of exactly identical windows of length -w <INT> that exist between two prealigned nucleotide sequences.

OPTIONS
  -d|--dir    [DIR] : dirname of directory containing alignments (fasta format)
  -o|--out   [FILE] : outfile prefix (default = mhom.table)
  -w|--window [INT] : window size range, comma delim, defined as: \"<SIZE>,<END>,<STEP>\" (default = 10,10,0)
  -h|--help         : this message

USAGE

\n";

my ($dirname, $outfile, $outprefix, $help);
my $windowrange = "10,10,1";

GetOptions (
  'd|dir=s'     => \$dirname,
  'out|o:s'     => \$outprefix,
  'window|w:s'  => \$windowrange,
  'help|h'      => \$help
);

die $usage if $help;
die $usage unless ($dirname);
if ($outprefix){
  $outfile = "$outprefix.mhom.table";
} else {
  $outfile = "mhom.table";
}
open (my $OUT, ">$outfile") or die "[ERROR] Cannot open $outfile: $!\n\n";
print STDERR "[INFO] Dirname: $dirname\n";
my @W = split (/\,/, $windowrange); ##get window sizes from string
print STDERR "[INFO] Window range: SIZE=$W[0],END=$W[1],STEP=$W[2]\n";

my @files = <$dirname/*fasta>;
print STDERR "[INFO] Number of files in $dirname: ".scalar(@files)."\n";

for ( my $window=$W[0];$window<=$W[1];$window+=$W[2] ) {
  print $OUT "$window";
  my $n=1;
  foreach my $file (nsort @files) {
    print STDERR "\r[INFO] Window size: $window; file \#$n";$| = 1;
    my $read_aln = Bio::AlignIO -> new(-file=>"$file", -format=>"fasta");
    my $aln = $read_aln->next_aln(); ##read aln
    my @seqs = $aln->each_seq(); ##get seqs from aln
    my @bits1 = ( $seqs[0]->seq() =~ /.{1,$window}/gs ); ##chop seq1 into bits1
    my @bits2 = ( $seqs[1]->seq() =~ /.{1,$window}/gs ); ##chop seq2 into bits2

    my $matches = 0;
    for my $i (0..$#bits1) {
      if ((length($bits1[$i])==$window) && (length($bits2[$i])==$window)) { ##only test if length(bit)==window
        $matches++ if $bits1[$i] eq $bits2[$i];
      }
    }
    print $OUT "\t$matches";
    $n++;
  }
  print $OUT "\n";
  print STDERR "\n";
}
close $OUT;
print "[INFO] Finished on ".`date`."\n";

__END__
