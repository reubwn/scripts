#!/usr/bin/env perl

## author: reubwn May 2019

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;

my $usage = "
SYNOPSIS:
  Strips trailing asterisks (stop codons; '*') from protein sequences
  Also removes sequences with internal stop codons

OPTIONS:
  -f|--fasta    [FILE] : fasta file [required]
  -i|--internal        : set option to also remove internal stop codons
  -h|--help            : prints this help message
\n";

my ($fasta,$internal,$help);

GetOptions (
  'f|fasta=s'     => \$fasta,
  'i|internal:s'  => \$internal,
  'h|help'        => \$help
);

die $usage if $help;
die $usage unless ($fasta);

my (%stripped,%removed);

my $in = Bio::SeqIO -> new ( -file => $fasta, -format => 'fasta' );
while (my $seq_obj = $in->next_seq()) {
  my $seq_string = $seq_obj->seq();
  if ($seq_string =~ m/\*$/) { ## has terminal stop codon
    $seq_string =~ s/\*$//; ## remove it
    if ($seq_string =~ m/\*/) { ## also has internal stop codon
      $removed{$seq_obj->display_id()}++; ## dont print it
    } else {
      print STDOUT ">" . $seq_obj->display_id() . "\n" . $seq_string . "\n";
      $stripped{$seq_obj->display_id()}++;
    }

  } elsif ($seq_string =~ m/\*/) { ## stop codon is internal; don't print
    $removed{$seq_obj->display_id()}++; ## dont print it
  } else {
    print STDOUT ">" . $seq_obj->display_id() . "\n" . $seq_string . "\n";
  }
}

print STDERR "[INFO] Number of seqs with trailing '*' stripped: " . (keys %stripped) . "\n";
print STDERR "[INFO] Number of seqs with internal '*' removed: " . (keys %removed) . "\n";
print STDERR "[INFO] Done " . `date`;
