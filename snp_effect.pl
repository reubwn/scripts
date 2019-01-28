#!/usr/bin/env perl

## Author: reubwn Jan 2019

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "
snp_effect.pl
=============

Calculates simple SNP effects and decomposes into synonymous/nonsynonymous,
noncoding, etc.

USAGE:
snp_effect.pl -i <vcf_file> -g <gff_file> -f <reference> [-o PREFIX]

OPTIONS:
  -i|--vcf   [FILE] : VCF file of variants [required]
  -g|--gff   [FILE] : GFF file of genes and/or other features [required]
  -f|--fasta [FILE] : fasta file of contigs [required]
  -o|--out   [STR]  : output filename [default: snp_effect.txt]
  -h|--help         : prints this help message
\n\n";

## params with defaults
my $prefix = "snp_effect.txt";

## other args
my ($vcf_file,$gff_file,$fasta,$help,$outfile);

GetOptions (
  'i|vcf=s'         => \$vcf_file,
  'f|fasta=s'       => \$fasta,
  'g|gff_file|g=s'  => \$gff_file,
  'o|out:s'         => \$prefix,
  'h|help'          => \$help,
);

die $usage if $help;
die $usage unless ($vcf_file && $gff_file && $fasta);
