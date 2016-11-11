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
  Takes a BLAST vs self file and returns a list of sequences with %-identity and
  %-query coverage over given thresholds. Use for finding transcripts that are
  very similar to other transcripts within a set (eg., duplicates or fragments).
  Will automatically discount Augustus-style alternative transcripts by looking
  for the .t1 and .t2 suffixes.

OUTPUTS:
  A list of the **shorter** of any highly similar sequences it finds.

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
  'blast|b=s'  => \$blast,
  'fastas|f=s' => \$fasta,
  'perc|p:f'   => \$perc,
  'qcov|q:f'   => \$qcov,
  'out|o:s'    => \$out,
  'help|h'     => \$help,
);

die $usage if $help;
die $usage unless ($fasta && $blast);

my $OUT;
if ($out) {
  open ($OUT, ">$out") or die "$!\n";
} else {
  $OUT = \*STDOUT;
}

## get seq lengths
my %seq_hash;
my $in = Bio::SeqIO->new ( -file => $fasta, -format => "fasta" );
while ( my $seq_obj = $in->next_seq() ){
  $seq_hash{($seq_obj->display_id())} = ($seq_obj->length());
}
print "Read in ".commify(scalar(keys %seq_hash))." sequences from $fasta \n\n";

## parse BLAST file.
## Assumes perc_identity is in $F[2] and qcovus is in $F[12] (correct if -outfmt '6 std qcovus')
my %filter_hash;
open (my $FILE, $blast) or die "$!\n";
while (<$FILE>) {
  my @F = split (/\t/, $_);
  ## account for possibility of alternative transcripts
  my @a = split (/\./, $F[0]);
  my @b = split (/\./, $F[1]);
  if ($a[0] eq $b[0]) {
    next; ## q and s are same gene
  } else {
    if (($F[2] >= $perc) and ($F[12] >= $qcov)) {
      my $qlen = $seq_hash{$F[0]};
      my $slen = $seq_hash{$F[1]};
      if ($slen <= $qlen) {
        $filter_hash{$F[1]}++;
      } else {
        $filter_hash{$F[0]}++;
      }
    }
  }
}
close $FILE;

## print results
foreach (nsort keys %filter_hash) {
  print $OUT "$_\n";
}
