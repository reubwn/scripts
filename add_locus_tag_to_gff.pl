#!/usr/bin/env perl

## author: reubwn Aug 2020

use strict;
use warnings;

use Getopt::Long;

my $usage = "
SYNOPSIS:
  Add protein_id and transcript_id descriptors to GFF file.
  The IDs are constructed using the locus tag and gene IDs, thus:
    protein_id = 'gnl|LOCUS_TAG|GENE_ID'
    transcript_id = 'gnl|LOCUS_TAG|GENE_ID_mrna'

  The attribute 'pseudo=true' can be applied to 'broken' genes, where the translation
  includes an internal stop codon probably due to assembly or annotation error. See
  https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/ for more info.

OPTIONS:
  -g|--gff         [FILE] : Input GFF file [required]
  -o|--out         [FILE] : Output filename [tagged.gff]
  -l|--locus_tag [STRING] : Given locus_tag for ID construct ['WGS:XXXX']
  -z|--zeroes       [INT] : Number of leading zeroes in numerator [6]
  -s|--suffix    [STRING] : Suffix for transcript ID ['']
  -p|--pseudo      [FILE] : Additionally annotate gene IDs in FILE with 'pseudo=true'
  -h|--help               : Prints this help message
\n";

my ($gff_file,$pseudo_file,$help);
my $gff_outfile = "tagged.gff";
my $leading_zeroes = 6;
my $locus_tag_prefix = "WGS:XXXX";
my $mrna_suffix = "";

GetOptions (
  'g|gff=s'        => \$gff_file,
  'o|outfile:s'    => \$gff_outfile,
  'l|locus_tag:s'  => \$locus_tag_prefix,
  'z|zeroes:i'     => \$leading_zeroes,
  's|suffix:s'     => \$mrna_suffix,
  'p|pseudo:s'     => \$pseudo_file,
  'h|help'         => \$help
);

die $usage if $help;
die $usage unless ($gff_file);

## locus tag enumerator
my $enumerator = 1;

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
  ## add locus_tag to gene features
  if ( $F[2] eq "gene" ) {
    my $locus_tag = $locus_tag_prefix . "_" . sprintf("%0".$leading_zeroes."d",$enumerator);
    print $OUT $line . ";" . "locus_tag=$locus_tag\n";
    $enumerator++;
  }
  ## add transcript_id to mRNA features
  elsif ($F[2] eq "mRNA") {
    my $ID = $1 if ($F[8] =~ m/ID=(.+?)(;|$)/); ## inherit ID from ID
    my $transcript_id = "transcript_id=gnl|$locus_tag_prefix|$ID" . $mrna_suffix; ## default is for transcript_id == protein_id
    print $OUT join (";", $line, $transcript_id) . "\n";
  }
  ## add protein_id to CDS features
  elsif ($F[2] eq "CDS") {
    my $ID = $1 if ($F[8] =~ m/Parent=(.+?)(;|$)/); ## inherit ID from Parent
    my $protein_id = "protein_id=gnl|$locus_tag_prefix|$ID";
    print $OUT join (";", $line, $protein_id) . "\n";
  } else {
    print $OUT "$line\n";
  }
}
close $IN;
close $OUT;

print STDERR "[INFO] Done " . `date`;

__END__
