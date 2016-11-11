#!/usr/bin/env perl

## author: reubwn October 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
SYNOPSIS:

OUTPUTS:

OPTIONS:
  -b|--blast [FILE]  : BLAST file created using \"-outfmt \'6 std qcovus\'\" [required]
  -f|--fasta [FILE]  : fasta file [required]
  -p|--perc  [FLOAT] : percentage identity threshold [default: 99\%]
  -q|--qcov  [FLOAT] : percentage query coverage threshold [default: 99\%]
  -o|--out   [STR]   : output filename [default: print to STDOUT]
  -h|--help          : prints this help message

EXAMPLES:

\n";

my ($blast,$fasta,$out,$help);
my $perc = 99;
my $qcov = 99;

GetOptions (
  'blast|b=s'    => \$blast,
  'fastas|f=s'    => \$fasta,
  'perc|p:f'    => \$perc,
  'qcov|q:f'      => \$qcov,
  'out|o:s'       => \$out,
  'help|h'        => \$help,
);

die $usage if $help;
die $usage unless ($fasta && $blast);
