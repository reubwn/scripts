#!/usr/bin/env perl

## author: reubwn May 2019

use strict;
use warnings;

use Getopt::Long;

my $usage = "
SYNOPSIS:
  Takes a list of UniProt accessions and returns gene name and taxonomy information.

OPTIONS:
  -i|--input [FILE] : query file of UniProt accessions (1 per line) [required]
  -p|--prot         : fetch protein fasta rather than text info
  -h|--help         : prints this help message
\n";

my ($infile,$prot,$help);

GetOptions (
  'i|input=s' => \$infile,
  'p|prot'    => \$prot,
  'h|help'    => \$help,
);

die $usage if $help;
#die $usage unless ($input);

my $IN; ## can read from STDIN
if ($infile){
  open ($IN, $infile) or die $!;
} else {
  $IN = *STDIN;
}

if ($prot) {
  while (my $string = <$IN>) {
    chomp $string;
    my @a = split (",", $string);
    foreach (@a) {
      print STDERR "[INFO] Fetching $_\n";
      `curl -sL "http://www.uniprot.org/uniprot/$_.fasta"`;
    }
  }
} else {
  while (my $string = <$IN>) {
    chomp $string;
    my @a = split (",", $string);
    foreach (@a) {
      `curl -sL 'http://www.uniprot.org/uniprot/$_.txt'`;
    }
  }
}

close $IN;
