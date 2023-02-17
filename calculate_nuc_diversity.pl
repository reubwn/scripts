#!/usr/bin/env perl

## reubwn Feb 23

use strict;
use warnings;
use Getopt::Long;

use File::Path qw( rmtree );
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

## RESULTS hash
my %RESULTS;
my %gene_lengths;

## parse populations files
open (my $POPS, $pops_file) or die $!;
my %pops;
while (my $pop = <$POPS>) {
  chomp $pop;
  $pops{$pop}++;
}
close $POPS;
print STDERR "[INFO] Number of populations in '$pops_file': ".scalar(keys %pops)."\n";

## make dir to store regions files
my $regions_dir = $outprefix."_regions";
if ( -d $regions_dir ) {
  print STDERR "[INFO] Emptying '$regions_dir'...\n";
  rmtree $regions_dir;
  mkdir $regions_dir;
} else {
  print STDERR "[INFO] Making '$regions_dir'\n";
  mkdir $regions_dir;
}

## open genes file
open (my $GENES, $genes_file) or die $!;

## grep gene IDs from GFF to generate regions.txt file
while (my $gene = <$GENES>) {
  chomp $gene;
  print STDERR "[INFO] > Gene: '$gene'\n";
  ## generate regions.txt file for BCF filtering
  open (my $REGIN, "grep $gene $gff_file |") or die $!;
  open (my $REGOUT, ">$regions_dir/$gene.regions.txt") or die $!;
  my $gene_length;
  while (my $line = <$REGIN>) {
    chomp $line;
    my @F = split (/\s+/, $line);
    if ($F[2] eq "CDS") {
      print $REGOUT join("\t", $F[0],$F[3],$F[4])."\n";
      $gene_length += (($F[4] - $F[3]) + 1); ## GFF coords are 1-based and inclusive
    }
  }
  close $REGIN;
  close $REGOUT;

  ## store total CDS length
  $gene_lengths{$gene} = $gene_length;

  ## iterate thru pops
  foreach my $pop (nsort keys %pops) {
    print STDERR "[INFO] -> Population: '$pop'\n";
    my $sum_pi = 0;

    ## get stats on SNPs in gene region, per population
    open (my $STATS, "bcftools view -R $regions_dir/$gene.regions.txt -S $samples_path/$pop.txt $vcf_file --min-ac=1 --no-update | bcftools stats - | grep \"\^SN\" |") or die $!;
    while (my $line = <$STATS>) {
      chomp $line;
      my @F = split (/\s+/, $line);
      $RESULTS{$gene}{$pop}{num_samples} = $F[-1] if $.==1;
      $RESULTS{$gene}{$pop}{num_snps} = $F[-1] if $.==4;
      $RESULTS{$gene}{$pop}{num_mnps} = $F[-1] if $.==5;
      $RESULTS{$gene}{$pop}{num_indels} = $F[-1] if $.==6;
      $RESULTS{$gene}{$pop}{num_multiallelic} = $F[-1] if $.==9;

    }
    close $STATS;

    ## skip if none
    if ( $RESULTS{$gene}{$pop}{num_snps} > 0 ) {
      print STDERR "[INFO] --> Found $RESULTS{$gene}{$pop}{num_snps} variant sites!\n";
      ## run vcftools --site-pi
      if ( system("bcftools view -R $regions_dir/$gene.regions.txt -S $samples_path/$pop.txt $vcf_file | vcftools --vcf - --out $gene.$pop --site-pi 2>/dev/null") != 0 ) {
        print STDERR "[INFO] --> Problem with vcftools command!\n";
      } else {
        print STDERR "[INFO] --> Ran vcftools successfully\n";
        open (my $SITESPI, "cut -f3 $gene.$pop.sites.pi |") or die $!;
        while (my $pi = <$SITESPI>) {
          next if $. == 1;
          chomp $pi;
          $sum_pi += $pi;
        }
        close $SITESPI;
        $RESULTS{$gene}{$pop}{pi} = ($sum_pi/$gene_lengths{$gene});
      }
    } else {
      print STDERR "[INFO] --> No variants found\n";
    }
  }
}
close $GENES;

foreach my $gene (nsort keys %RESULTS) {
  my %pops = %{$RESULTS{$gene}};
  foreach my $pop (nsort keys %pops) {
    print STDOUT join ("\t", $gene, $gene_lengths{$gene}, $pop, $pops{num_samples}, $pops{num_snps}, $pops{num_mnps}, $pops{num_indels}, $pops{num_multiallelic}, $pops{pi}) . "\n";
  }
}
