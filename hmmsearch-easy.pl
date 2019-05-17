#!/usr/bin/env perl

## author: reubwn May 2019

use strict;
use warnings;

use Getopt::Long;

my $usage = "
SYNOPSIS:
  Runs hmmsearch vs fasta database, and returns alignment (clustal) including
  seeds in the HMM file.

OPTIONS:
  -q|--query  [FILE]   : query file (HMM format) [required]
  -d|--db     [FILE]   : sequence database (fasta format) [required]
  -a|--ali    [FILE]   : query file (aligned .sto or .fasta) [required]
  -e|--evalue [STRING] : hmmsearch reporting evalue threshold [1e-5]
  -p|--prefix [STRING] : output file prefix ['hmmsearch-easy']
  -h|--help            : prints this help message
\n";

my ($query,$db,$ali,$help);
my $evalue = "1e-5";
my $prefix = "hmmsearch-easy";

GetOptions (
  'q|query=s'    => \$query,
  'd|db=s'    => \$db,
  'a|ali=s' => \$ali,
  'e|evalue:s'      => \$evalue,
  'p|prefix:s' => \$prefix,
  'help|h'        => \$help,
);

die $usage if $help;
die $usage unless ($query && $db && $ali);

if (system("hmmsearch -h &> /dev/null")!=0) {
  die "[ERROR] Is hmmsearch installed and in \$PATH?\n";
}

if (system("fastaqual_select.pl &>/dev/null")!=0) {
  die "[ERROR] Is fastaqual_select.pl installed and in \$PATH?\n";
}

print STDERR "[INFO] Query HMM: $query\n";
print STDERR "[INFO] Query alignment: $ali\n";
print STDERR "[INFO] Database: $db\n";

## run hmmsearch
`hmmsearch -o $prefix.hmmsearch.$evalue.out --tblout $prefix.hmmsearch.$evalue.tblout --noali -E $evalue $query $db`;

## get hits to faa
`awk '{print \$1}' $prefix.hmmsearch.$evalue.tblout | fastaqual_select.pl -f $db -i - > $prefix.hmmsearch.$evalue.out.faa`;

## run hmmalign
`hmmalign -o $prefix.hmmalign.sto --mapali $ali --trim $query $prefix.hmmsearch.$evalue.out.faa`;

## convert aligned *.sto to clustal
`esl-reformat clustal $prefix.hmmalign.sto > $prefix.hmmalign.clustal`;

print STDERR "[INFO] Done " . `date`;
