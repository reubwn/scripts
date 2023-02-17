#!/usr/bin/env perl

## reubwn Feb 23

use strict;
use warnings;
use Getopt::Long;

use File::Basename;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Calculate per-gene per-population avg. nucleotide diversity (pi) from GFF and VCF.

OPTIONS [*required]
  -v|--vcf      *[FILE] : VCF/BCF file
  -g|--genes    *[FILE] : TXT list of gene IDs to be analysed (linking to the 'mRNA' field of the GFF)
  -G|--gff      *[FILE] : GFF annotation file
  -p|--pops     *[FILE] : TXT list of population/site IDs to subset VCF
  -s|--samples  *[PATH] : path/to/dir containing samples files that link sample ID to population grouping
  -o|--out       [STR]  : outfiles prefix ('pi')
  -h|--help             : print this message
\n";

my ($vcf_file, $genes_file, $gff_file, $pops_file, $samples_path, $help);
my $outprefix = "pi";

GetOptions (
  'v|vcf=s'     => \$vcf_file,
  'g|genes=s'   => \$genes_file,
  'G|gff=s'     => \$gff_file,
  'p|pops=s'    => \$pops_file,
  's|samples=s' => \$samples_path,
  'o|out:s'     => \$outprefix,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $vcf_file && $genes_file && $pops_file && $samples_path );

## open genes and populations files
open (my $GENES, $genes_file) or die $!;
open (my $POPS, $pops_file) or die $!;

## grep gene IDs from GFF to generate regions.txt file
while (my $gene = <$GENES>) {
  chomp $gene;
  print STDERR "[INFO] Gene: '$gene'\n";
  ## generate regions.txt file for BCF filtering
  open (my $REGIN, "grep $gene $gff_file |") or die $!;
  open (my $REGOUT, ">$gene.regions.txt") or die $!;
  while (my $line = <$REGIN>) {
    chomp $line;
    my @F = split (/\s+/, $line);
    if($F[2] eq "CDS"){
      print $REGOUT join("\t", $F[0],$F[3],$F[4]) . "\n";
    }
  }
  close $REGIN;
  close $REGOUT;
  # `grep $gene $gff_file | perl -lane 'if($F[2]eq"CDS"){print join("\t",$F[0],$F[3],$F[4])}' > ${gene}.regions.txt`;

  ## iterate thru pops
  while (my $pop = <$POPS>) {
    chomp $pop;
    print STDERR "[INFO] Population: '$pop'\n";
    my $num_variants = `bcftools view -R $gene.regions.txt -S $samples_path/$pop.txt $vcf_file | grep -v "^#" | wc -l`;
    ## check if there are any variants intersecting with gene/pop
    if ( $num_variants > 0 ) {
      ## run vcftools --site-pi
      print STDERR "[INFO] Found $num_variants variant lines for '$gene' in '$pop'\n";
      if ( system("bcftools view -R $gene.regions.txt -S $samples_path/$pop.txt $vcf_file | vcftools --vcf - --out $gene.$pop --site-pi") != 0 ) {
        print STDERR "[INFO] Problem with vcftools command!\n";
      } else {
        print STDERR "[INFO] Ran vcftools successfully\n";
      }
    } else {
      print STDERR "[INFO] No variants found for '$gene' in '$pop'\n";
    }
  }
  close $POPS;

}
close $GENES;
