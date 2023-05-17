#!/usr/bin/env perl

## reubwn May 23

use strict;
use warnings;
use Getopt::Long;

use File::Path qw( rmtree );
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Convert OrthoFinder Orthogroups.txt file to phyletic presence/absence matrix.

OPTIONS [*required]
  -i|--in      *[FILE] : Orthogroups.txt file
  -p|--path    *[PATH] : path to protein files
  -e|--ext     *[STR]  : protein file extension to glob ("*.faa")
  -o|--out      [STR]  : outfiles prefix ('phyletic')
  -h|--help            : print this message
\n";

my ($orthogroups_file, $proteins_path, $help);
my $extenstion = "faa";
my $outprefix = "phyletic";

GetOptions (
  'v|vcf=s'     => \$vcf_file,
  'g|genes=s'   => \$genes_file,
  'G|gff=s'     => \$gff_file,
  'p|pops=s'    => \$pops_file,
  's|samples=s' => \$samples_path,
  'o|out:s'     => \$outprefix,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $vcf_file && $genes_file && $pops_file && $samples_path );
