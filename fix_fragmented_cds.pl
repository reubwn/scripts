#!/usr/bin/env perl

## reubwn December 2022

use strict;
use warnings;
use Getopt::Long;

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Fix fragmented CDS from poor quality genome annotations using inference from orthology.

OPTIONS [*required]
  -a|--aa          *[FILE] : target aa sequences (fasta)
  -d|--db          *[DIR]  : dir of database proteome files (fasta)
  -g|--orthogroups *[FILE] : OrthoGroups.txt file (from OrthoFinder)
  -n|--dna          [FILE] : target DNA sequences (fasta)
  -o|--out          [STR]  : outfile suffix
  -l|--logfile             : print stats to logfile [no]
  -h|--help                : print this message
\n";

my ($aa_file, $dna_file, $db_dir, $orthogroups_file, $logfile, $help);
my $outsuffix = "trimmed";

GetOptions (
  'a|aa=s'      => \$aa_file,
  'd|db=s'      => \$db_dir,
  'g|orthogroups=s' => \$orthogroups_file,
  'n|dna:s'     => \$dna_file,
  'o|out:s'     => \$outsuffix,
  'l|logfile'   => \$logfile,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $aa_file && $orthogroups_file && $db_dir );
