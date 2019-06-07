#!/usr/bin/env perl

## author: reubwn May 2019

use strict;
use warnings;

use Getopt::Long;

my $usage = "
SYNOPSIS:
  Quick download GenBank from accession.

OPTIONS:
  -a|--acc   [STRING] : GenBank accession number [required]
  -f|--fasta          : also pull region as nucleotide fasta
  -p|--prot           : also pull the CDS as proteins
  -c|--cds            : also pull the CDS as nucleotides
  -h|--help           : prints this
\n";

my ($accession_string,$fasta,$prot,$cds,$help);

GetOptions (
  'a|acc=s' => \$accession_string,
  'f|fasta' => \$fasta,
  'p|prot'  => \$prot,
  'c|cds'   => \$cds,
  'h|help'  => \$help,
);

die $usage if $help;
die $usage unless ($accession_string);

my @array = split (/,/, $accession_string);
foreach my $accession (@array) {
  chomp $accession;
  print "[INFO] Accession number is: $accession\n";
  if ( system("curl -s  'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=gb&retmode=txt'>$accession.gbk") != 0 ) {
    die "[ERROR] Problem with eutils command for GenBank!\n";
  }

  if ($fasta) {
    if ( system("curl -s  'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta&retmode=txt'>$accession.fasta") != 0 ) {
      die "[ERROR] Problem with eutils command for fasta!\n";
    }
  }

  if ($prot) {
    if ( system("curl -s  'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta_cds_aa&retmode=txt'>$accession.CDS.faa") != 0 ) {
      die "[ERROR] Problem with eutils command for CDS faa!\n";
    }
  }

  if ($cds) {
    if ( system("curl -s  'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta_cds_na&retmode=txt'>$accession.CDS.fna") != 0 ) {
      die "[ERROR] Problem with eutils command for CDS fna!\n";
    }
  }
}

print STDERR "[INFO] Done " . `date`;
