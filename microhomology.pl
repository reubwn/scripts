#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
SYNOPSIS

OPTIONS
  -i|--in        [FILE]: infile
  -o|--out       [FILE]: outfile (default = STDOUT)
  -w|--window    [INT] : window size (default = 10)
  -h|--help            : this message

USAGE

\n";

my ($infile, $outfile, $help);
my $window = 10;

GetOptions (
  'in|i=s'      => \$infile,
  'out|o:s'     => \$outfile,
  'window|w:i'  => \$window,
  'help|h'      => \$help
);

die $usage if $help;
#die $usage unless $infile;
print STDERR "[INFO] File: $infile\n";
print STDERR "[INFO] Window size: $window\n";

my $seq1 = "ATCGATCGATCGATCGATCG";
my $seq2 = "AACGA-CGATCGANCGATCG";
my $matches;
my (%substr1,%substr2);

my @bits1 = ( $seq1 =~ /.{1,$window}/gs );
my @bits2 = ( $seq2 =~ /.{1,$window}/gs );
print "@bits1\n";

for my $i (0..$#bits1) {
  if ((length($bits1[$i])==$window) && (length($bits2[$i])==$window)) {
    print "1:$bits1[$i]\n2:$bits2[$i]\n";
    print "Match!\n" if $bits1[$i] eq $bits2[$i];
    $matches++ if $bits1[$i] eq $bits2[$i];
  } else {print "Window too short!\n"}
}
print "Total matches: $matches\n";
