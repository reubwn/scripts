#!/usr/bin/env perl

## author: reubwn June 2017

use strict;
use warnings;
use Getopt::Long;
use List::Util qw /sum/;
use Sort::Naturally;

my $usage = "
SYNOPSIS

OUTPUT

OPTIONS
  -a|--bam [FILE] : BAM file
  -b|--bed        : windows in BED format
  -g|--genome     : genome file that defines the expected chromosome order in the input files
  -o|--out        : outfile name
  -h|--help       : print this message

USAGE

\n";

my ($bam,$bed,$genome,$help,$debug);
my $n = 1;

GetOptions (
  'a|bam=s'  => \$bam,
  'b|bed=s'  => \$bed,
  'g|genome=s' => \$genome,
  'h|help'   => \$help,
  'd|debug'  => \$debug
);

die $usage if $help;
die $usage unless ($bam && $bed && $genome);

print STDERR "[INFO] BAM file: $bam\n";
print STDERR "[INFO] Windows file: $bed\n";
print STDERR "[INFO] Genome file: $genome\n";

open (my $BED, $bed) or die "[ERROR] Cannot open $bed: $!\n";
while (my $window = <$BED>) {
  #print STDERR "\r[INFO] Working on window \#$n: $_";$| = 1;
  `printf "$window" > tmp.bed`;
  my ($SAM,$total,$same,$insert,@insert_arr,$insert_avg);
  open($SAM, "bedtools intersect -sorted -g $genome -a $bam -b tmp.bed | samtools view - |");#`bedtools intersect -sorted -g $genome -a $bam -b <(printf "$_") | samtools view - | perl -lane 'if($F[6]eq"="){if($F[8]>500){$insert++};$same++;$total++}else{$total++}END{print "$total\t$same\t".($same/$total)."\t$insert\t".($insert/$total)}'`;
  while (<$SAM>) {
    my @F = split (/\s+/, $_);
    if ($F[6] eq "=") { ##mate on same scaffold
      $same++;
      $total++;
      $insert++ if ($F[8] > 500);
      push (@insert_arr, $F[8]);
    } else {
      $total++;
    }
  }
  close $SAM;
  chomp($window);
  print STDOUT join (
    "\t",
    $window,
    $total,
    $same,
    ($same/$total),
    $insert,
    ($insert/$total),
    (sum(@insert_arr)/scalar(@insert_arr)),
    "\n"
  );
}

close $BED;
print STDERR "\n[INFO] Finished on ".`date`."\n";
