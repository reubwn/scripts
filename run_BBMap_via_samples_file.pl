#!/usr/bin/env perl

## author: reubwn Sep 2020

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;
use Sort::Naturally;

my $usage = "
SYNOPSIS:
  Run BBMap from Trinity samples file.

OPTIONS:
  -g|--genome       [FILE] : Target genome file [required]
  -s|--samples_file [FILE] : Trinity samples file [required]
  -t|--threads       [INT] : Number of threads [4]
  -m|--memory        [INT] : Max memory in Gb [62]
  -h|--help                : Prints this help message
\n";

my ($genome_file,$samples_file,$help);
my $threads = 4;
my $memory = 62;

GetOptions (
  'g|genome=s'       => \$genome_file,
  's|samples_file=s' => \$samples_file,
  't|threads:i'      => \$threads,
  'm|memory:i'       => \$memory,
  'h|help'           => \$help
);

die $usage if $help;
die $usage unless ($genome_file && $samples_file);

## parse samples file
my %samples_hash;
open (my $SAMPLES, $samples_file) or die $!;
while (my $line=<$SAMPLES>) {
  chomp $line;
  my @F = split (m/\s+/, $line);
  $samples_hash{$F[1]}{R1} = $F[2];
  $samples_hash{$F[1]}{R2} = $F[3];
}
close $SAMPLES;

# print Dumper \%samples_hash;

## run bbmap sequentially
foreach my $sample (nsort keys %samples_hash) {
  print STDERR "[INFO] Running BBMap on $sample...\n";
  `bbmap.sh -Xmx${memory}g nodisk=t ref=$genome_file in1=$samples_hash{$sample}{R1} in2=$samples_hash{$sample}{R2} maxindel=50k threads=$threads outm=${sample}_mapped.sam.gz statsfile=${sample}_mapped.sam.stats`;
  print STDERR "[INFO] Done " . `date`;
}

print STDERR "[INFO] All done! " . `date`;
