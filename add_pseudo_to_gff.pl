#!/usr/bin/env perl

## author: reubwn Aug 2020

use strict;
use warnings;

use Getopt::Long;

my $usage = "
SYNOPSIS:
  The attribute 'pseudo=true' can be applied to 'broken' genes, where the translation
  includes an internal stop codon probably due to assembly or annotation error. See
  https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/ for more info.

OPTIONS:
  -g|--gff         [FILE] : Input GFF file [required]
  -o|--out         [FILE] : Output filename [tagged.gff]
  -p|--pseudo      [FILE] : Additionally annotate gene IDs in FILE with 'pseudo=true'
  -h|--help               : Prints this help message
\n";

my ($gff_file,$pseudo_file,$help,$debug);
my $gff_outfile = "pseudo_true.gff";

GetOptions (
  'g|gff=s'     => \$gff_file,
  'o|outfile:s' => \$gff_outfile,
  'p|pseudo:s'  => \$pseudo_file,
  'h|help'      => \$help,
  'd|debug'     => \$debug
);

die $usage if $help;
die $usage unless ($gff_file);

my %broken_genes;
if ( $pseudo_file ) {
  open (my $P, $pseudo_file) or die $!;
  while (<$P>) {
    chomp;
    $broken_genes{$_} = ();
  }
  close $P;
  print STDERR "[INFO] Number of genes to be annotated 'pseudo=true': " . scalar(keys %broken_genes) . "\n";
  print STDERR join (" ", keys %broken_genes) . "\n" if ( $debug );
}

## IO
open (my $IN, $gff_file) or die $!;
open (my $OUT, ">$gff_outfile") or die $!;
print STDERR "[INFO] Writing to '$gff_outfile'\n";

while (my $line = <$IN>) {
  chomp $line;

  if ($line =~ m/^#/) {
    print $OUT "$line\n";
    next;
  }
  $line =~ s/;$//; ## trim trailing ';' if present
  my @F = split (/\s+/, $line);

  ## find broken genes from list
  if ( $F[2] eq "gene" ) {
    my $ID = $1 if ($F[8] =~ m/locus_tag=(.+?)(;|$)/); ## inherit ID from locus_tag
    print STDERR "Found locus_tag: $ID\r" if ( $debug );
    if ( $broken_genes{'H9Q69_000736'} ) {
      ## gene is broken
      print STDERR "Found ID: $ID\n" if ( $debug );
      print $OUT $line . ";" . "pseudo=true\n"; ## only required for 'gene' feature
    } else {
      print $OUT "$line\n";
    }
  } else {
    print $OUT "$line\n";
  }
}
close $IN;
close $OUT;

print STDERR "[INFO] Done " . `date`;

__END__
