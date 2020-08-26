#!/usr/bin/env perl

## author: reubwn Aug 2020

use strict;
use warnings;

use File::Copy;
use Getopt::Long;

my $usage = "
SYNOPSIS:
  Add locus_tag prefix to .tbl file (output from GAG).
  Prints to STDOUT

OPTIONS:
  -t|--tbl         [FILE] : Input .tbl file (e.g. from GAG) [required]
  -l|--locus_tag [STRING] : Replace current locus_tag text with STRING ['ABCD']
  -z|--zeroes       [INT] : Number of leading zeroes in numerator [6]
  -b|--bak                : Make a backup (.bak file)
  -h|--help               : Prints this help message
\n";

my ($tbl_file,$bak,$help);
my $locus_tag_prefix = "ABCD";
my $leading_zeroes = 6;

GetOptions (
  't|tbl=s'        => \$tbl_file,
  'l|locus_tag:s'  => \$locus_tag_prefix,
  'z|zeroes:i'     => \$leading_zeroes,
  'b|bak'          => \$bak,
  'h|help'         => \$help
);

die $usage if $help;
die $usage unless ($tbl_file);

## create .bak
if ($bak) {
  `cp $tbl_file $tbl_file.bak`;
}

## replace locus_tags
my $n = 1;
open (my $IN, $tbl_file) or die $!;
open (my $OUT, ">$tbl_file.tmp") or die $!;
while (<$IN>) {
  chomp;
  if (m/locus_tag/) {
    print $OUT "\t\t\tlocus_tag\t$locus_tag_prefix\_".sprintf("%0".$leading_zeroes."d",$n)."\n";
    $n++;
  } else {
    print $OUT "$_\n";
  }
}
close $IN;
close $OUT;

## clean up
if ($bak) {
  `mv $tbl_file.tmp $tbl_file`;
  `gzip $tbl_file.bak`;
} else {
  `mv $tbl_file.tmp $tbl_file`;
}

print STDERR "Done. ".`date`;

__END__
