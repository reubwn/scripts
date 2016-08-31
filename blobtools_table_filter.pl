#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;

my $usage = "
Blobtools Table Filter
Parses blobtools table from \`blobtools view -r all -b\` and returns some statistics.

USAGE: btf.pl

OPTIONS:
  -i|--in       : blobtools *table.txt file [required]
  -s|--select   : primary search query... i.e., what are you looking for? [required]
                  (regex style case insensitive, Proteobacteria == proteobac etc...)
  -r|--rank     : taxonomic rank to at which to search query given by -s [default: P]
                  ('Su'=superkingdom, 'P'=phylum, 'O'=order, 'C'=class, 'F'=family, 'G'=genus, 'Sp'=species)
  -o|--outtype  : what to output [default: contig names]
                  ('N'=contig names, 'S'=sequences in fasta format [requires -f <fasta_file>], 'F'=full table entry)
  -g|--gc       : minimum GC proportion filter [default: 0]
  -a|--at       : minimum AT proportion filter [default: 0]
  -l|--length   : minimum length filter [default: none]
  -c|--coverage : minimum coverage filter [default: none]
  -f|--fasta    : sequences in fasta format
  -h|--help     : prints this help message

EXAMPLES:
  (1) Get all contig names for contigs labelled as 'Proteobacteria' at the phylum level, with G+C > 50\%:
      >> btf.pl -i blobtools.table.txt -s \"proteobac\" --gc 0.5 > proteobacterial_contigs.list
  (2) Get sequences of contigs labelled 'Trichoplax adhaerens' greater than 500 nt lenth:
      >> btf.pl -i blobtools.table.txt -s 'adhaerens' -r 'Sp' -o 'S' -l 500 > Tadhaerens_contigs.fasta
  (3) Get probably the mitochondrion sequences:
      >> btf.pl -i blobtools.table.txt --at 0.45 --cov 200
\n";

## params with defaults
my $rank = "P";
my $outtype = "N";

## other args
my ($in_file,$select,$gc,$at,$length,$coverage,$fasta,$help);

GetOptions (
'in|i=s'      => \$in_file,
'select|s=s'  => \$select,
'rank|r:s'    => \$rank,
'outtype|o:s' => \$outtype,
'gc|g:f'      => \$gc,
'length|l:i'  => \$length,
'coverage|c:f'=> \$coverage,
'fasta|f:s'   => \$fasta,
'rank|r:s'    => \$rank,
'help|h'      => \$help,
);

die $usage if $help;
die $usage unless $in_file;
