#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;
use Getopt::Long;
use List::Util qw /sum/;
use Sort::Naturally;

my $usage = "
SYNOPSIS

OUTPUT

OPTIONS
  -i|--collin [FILE] : collinearity file from MCScanX (annotated with Ks)
  -k|--ks            : only examine blocks with Ks <= this threshold
  -h|--help          : print this message

USAGE

\n";

my ($collinearity,$help,$debug);
my $ks = 0.5;

GetOptions (
  'collinearity|i=s' => \$collinearity,
  'kaks|k'           => \$ks,
  'help|h'           => \$help,
  'debug|d'          => \$debug
);

die $usage if $help;
die $usage unless ($collinearity);

print STDERR "[INFO] Collinearity file: $collinearity\n";
print STDERR "[INFO] Ks threshold: $ks\n";

my (%blocks1,%blocks2,%chroms1,%chroms2);

open (my $IN, $collinearity) or die "$!\n";
while (<$IN>) {
  next if $. == 1; ##skip header on line1
  my @F = split (/\s+/, $_);
  if ($F[11] <= $ks) {
    $chroms1{$F[0]} = $F[1]; ##key= blockname; val=chrom where block exists
    push ( @{$blocks1{$F[1]}}, $F[0] ); ##key=chrom1 name; val= @{ blocks on chrom1 name }
    $chroms2{$F[0]} = $F[2]; ##key= blockname; val=chrom where block exists
    push ( @{$blocks2{$F[2]}}, $F[0] ); ##key=chrom2 name; val= @{ blocks on chrom2 name }
  }
}
close $IN;

# foreach (nsort keys %blocks1) {
#   my @a = @{$blocks1{$_}};
#   print "$_\t@a\n";
# }

foreach my $block (nsort keys %chroms1) {
  #print "[INFO] Block #$block: $chroms1{$block}\n"; ##name of chrom where $block resides
  my @blocksOnSameChrom = @{ $blocks1{$chroms1{$block}} }; ##all blocks residing on the same chrom
  my @blocksOnHomologousChrom = @{ $blocks2{$chroms2{$block}} }; ##all other blocks residing on (what should be) the HOMOLOGOUS chrom
  my( $index1 )= grep { $blocksOnSameChrom[$_] == $block } 0..$#blocksOnSameChrom; ##get index of block in series of blocks on same chrom
  my( $index2 )= grep { $blocksOnHomologousChrom[$_] == $block } 0..$#blocksOnHomologousChrom; ##get index of HOMOLOGOUS block on HOMOLOGOUS chrom
  print "$block\t$chroms1{$block}\t$index1\t$chroms2{$block}\t$index2\t@blocksOnHomologousChrom\n";

}
