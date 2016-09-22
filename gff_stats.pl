#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
GFF_stats.pl
============

Print some basic metrics from a GFF file.

OPTIONS:
  -i|--in   : input GFF file
  -o|--out  : outfile prefix (default = STDOUT)
  -h|--help : prints this message
";

my ($in, $help);
my $out = STDOUT;

GetOptions (
  'in|i=s'                 => \$gff_file,
  'out|o:s'                => \$prefix,
  'help|h'                 => \$help,
);

die if $help;
die unless $in;

my ($transcript_name)

open (my $IN, "<", $in) or die $!;
while (<$IN>){
  next if /^\#/;
  my @F = split /\t/, $_;
  $transcript_name = $F[7] if $F[2] eq /transcript/;
  print $out $transcript_name."\n";
}
