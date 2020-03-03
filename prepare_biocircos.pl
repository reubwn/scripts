#!/usr/bin/env perl

## author: reubwn October 2019

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
Prints genome object in format for BioCircos input.

OPTIONS
  -f|--fasta  [FILE] : input fasta file [required]
  -n|--name    [STR] : genome 'name' ('myGenome')
  -l|--length  [INT] : length threshold (0)
  -h|--help          : prints this help message
\n";

my ($fasta,$length,$help);
my $name = "myGenome";

GetOptions (
  'f|fasta=s'  => \$fasta,
  'n|name:s'   => \$name,
  'l|length:i' => \$length,
  'h|help'     => \$help,
);

die $usage if $help;
die $usage unless $fasta;

my %seq_hash;

print "$name = list(\n";

my $in = Bio::SeqIO -> new ( -file => $fasta, -format => 'fasta' );
while (my $seq_obj = $in -> next_seq()) {
  next if $seq_obj->length() < $length;
  print "  '" . $seq_obj->display_id() . "' = " . $seq_obj->length();
  print ",\n" unless eof;
}

print "\n)\n";

__END__
