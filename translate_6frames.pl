#!/usr/bin/env perl

## Author: reubwn Jan 2019

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqUtils;
use Getopt::Long;

my $usage = "
translate_6frames.pl
=============

Translates a genome fasta into all 6 forward and reverse coding frames, e.g. for
a 'BLASTX-like' analysis using HMMER.

USAGE:
translate_6frames.pl -f <fasta> > FILE

OPTIONS:
  -f|--fasta [FILE] : fasta file of contigs [required]
  -h|--help         : prints this help message
\n\n";

## other args
my ($fasta,$help);

GetOptions (
  'f|fasta=s'       => \$fasta,
  'h|help'          => \$help,
);

die $usage if $help;
die $usage unless ($fasta);

## open fasta
my $in = Bio::SeqIO -> new ( -file => $fasta, -format => "fasta" );
my $stream = Bio::SeqIO->newFh( -fh => \*STDOUT, -format => "fasta" );

while ( my $seq_obj = $in->next_seq() ) {
  my @translated = Bio::SeqUtils->translate_6frames($seq_obj);
  foreach (@translated) {
    print $stream $_;
  }
}

print STDERR `date`;
