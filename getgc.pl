#!/usr/bin/env perl

## author: reubwn May 2019

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;

my $usage = "
SYNOPSIS:
  Count GC per input sequence.
  Prints to STDOUT

OPTIONS:
  -f|--fasta [FILE]   : fasta file [required]
  -t|--text  [STRING] : add TEXT as extra columns in output
  -h|--help           : prints this help message
\n";

my ($fasta,$text,$help);

GetOptions (
  'f|fasta=s' => \$fasta,
  't|text:s'  => \$text,
  'h|help'    => \$help,
);

die $usage if $help;
die $usage unless ($fasta);

my $in = Bio::SeqIO -> new ( -file => $fasta, -format => 'fasta' );
while ( my $seq_obj = $in -> next_seq() ) {
  ## count G+C (inc soft masked gc)
  my $gc_count = $seq_obj->seq() =~ tr/GCgc//;
  if ($text) {
    print STDOUT join ("\t", $seq_obj->display_id, ($gc_count/$seq_obj->length()), $text ) . "\n";
  } else {
    print STDOUT join ("\t", $seq_obj->display_id, ($gc_count/$seq_obj->length()) ) . "\n";
  }
}

__END__
