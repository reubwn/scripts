#!/usr/bin/env perl

## reubwn Feb 23

use strict;
use warnings;
use Getopt::Long;

use File::Basename;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Calculate per-gene per-population avg. nucleotide diversity (pi) from GFF and VCF.

OPTIONS [*required]
  -v|--vcf      *[FILE] : VCF/BCF file
  -g|--genes    *[FILE] : list of gene IDs to be analysed (linking to the 'mRNA' field of the GFF)
  -G|--gff      *[FILE] : GFF annotation file
  -s|--samples  *[FILE] : TAB delim samples file linking sample ID to population grouping (sampleID,popID)
  -o|--out       [STR]  : outfile prefix ('nuc_diversity.txt')
  -h|--help             : print this message
\n";

my ($vcf_file, $genes_file, $gff_file, $samples_file, $help);
my $outprefix = "nuc_diversity";

GetOptions (
  'v|vcf=s'     => \$vcf_file,
  'g|genes=s'   => \$genes_file,
  'G|gff=s'     => \$gff_file,
  's|samples=s' => \$samples_file,
  'o|out:s'     => \$outprefix,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $vcf_file && $genes_file && $samples_file );
