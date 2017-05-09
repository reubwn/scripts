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
  Uses a sliding window incrementing by 1, and scales number of identical windows by the length of each alignment.

OPTIONS
  -d|--dir    [DIR] : dirname of directory containing alignments (fasta format)
  -o|--out   [FILE] : outfile prefix (default = mhom.table)
  -w|--window [INT] : window size range, comma delim, defined as: \"<SIZE>,<END>,<STEP>\" (default = 10,10,1)
  -h|--help         : this message

USAGE
  (1) To calculate number of idenical windows from window size 1 to 40, increasing window size by 2 each round:
  >> microhomology.pl -d <alignments/> -w \"1,40,2\" -o out
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
print STDERR "[INFO] Window range: MIN=$W[0],MAX=$W[1],STEP=$W[2]\n";

my @files = <$dirname/*>;
my @filenames = map {$_ =~ s/.+\///r} @files;
print $OUT join("\t", "WINDOW", (nsort @filenames), "\n");
print STDERR "[INFO] Number of files in $dirname: ".scalar(@files)."\n";

## iterate across window range; from MIN to MAX testing every STEP
for ( my $window=$W[0];$window<=$W[1];$window+=$W[2] ) {
  print $OUT "$window";
  my $n=1;
  foreach my $file (nsort @files) {
    print STDERR "\r[INFO] Window size: $window; file \#$n";$| = 1;
    my $read_aln = Bio::AlignIO -> new(-file=>"$file", -format=>"fasta");
    my $aln = $read_aln->next_aln(); ##read aln
    my $aln_length = $aln->length(); ##aln length for scaling
    my @seqs = $aln->each_seq(); ##get seqs from aln; this does print the gaps too
    my $matches = 0;

    ## implement sliding window, always incrementing by 1:
    ## number of windows given sliding window +1 is given by: (aln length)-(window size)+1
    for ( my $j=1;$j<=($aln_length-$window+1);$j++ ) {
      ## count if subseqs from start ($j) to end ($j+$window-1) are identical:
      if ( ($seqs[0]->subseq($j,($j+$window-1))) eq ($seqs[1]->subseq($j,($j+$window-1)))) {
        $matches++;
      }
    }
    ## below code is how to chop sequences into discrete windows; obsoleted
    # my @bits1 = ( $seqs[0]->seq() =~ /.{1,$window}/gs ); ##chop seq1 into bits1
    # my @bits2 = ( $seqs[1]->seq() =~ /.{1,$window}/gs ); ##chop seq2 into bits2
    #
    # my $matches = 0;
    # for my $i (0..$#bits1) {
    #   if ((length($bits1[$i])==$window) && (length($bits2[$i])==$window)) { ##only test if length(bit)==window
    #     $matches++ if $bits1[$i] eq $bits2[$i];
    #   }
    # }

    ## print the number of matches / alignment length, as a scaling factor to account for length differences in genes:
    print $OUT "\t".($matches/$aln_length);
    $n++;
  }
  print $OUT "\n";
  print STDERR "\n";
}
close $OUT;
print "[INFO] Finished on ".`date`."\n";

__END__
