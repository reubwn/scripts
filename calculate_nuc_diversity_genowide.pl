#!/usr/bin/env perl

## reubwn July 23

use strict;
use warnings;
use Getopt::Long;

use File::Path qw( rmtree );
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Calculate genome wide nucleotide diversity (pi) across provided samples.

OPTIONS [*required]
  -i|--in      **[FILE] : FOFN list of VCF/BCF files to analyse, or...
  -v|--vcf     **[FILE] : single VCF/BCF file
  -b|--bed      *[FILE] : BED file of regions to be included in analysis
  -p|--pops     *[FILE] : TXT list of population/site IDs to subset VCF
  -s|--samples  *[PATH] : path/to/dir containing samples files that link sample ID to population grouping
  -o|--out       [STR]  : outfiles prefix ('result')
  -h|--help             : print this message
\n";

my ($vcf_fofn, $vcf_file, $bed_file, $pops_file, $samples_path, $help);
my $outprefix = "result";

GetOptions (
  'i|in:s'      => \$vcf_fofn,
  'v|vcf:s'     => \$vcf_file,
  'b|bed=s'     => \$bed_file,
  'p|pops=s'    => \$pops_file,
  's|samples=s' => \$samples_path,
  'o|out:s'     => \$outprefix,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $vcf_fofn || $vcf_file );
die $usage unless ( $pops_file && $samples_path );

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

## make dir to store intermediate files
my $regions_dir = $outprefix."_regions";
if ( -d $regions_dir ) {
  print STDERR "[INFO] Emptying '$regions_dir'...\n";
  rmtree $regions_dir;
  mkdir $regions_dir;
} else {
  print STDERR "[INFO] Making '$regions_dir'\n";
  mkdir $regions_dir;
}

my $pi_dir = $outprefix."_sites_pi";
if ( -d $pi_dir ) {
  print STDERR "[INFO] Emptying '$pi_dir'...\n";
  rmtree $pi_dir;
  mkdir $pi_dir;
} else {
  print STDERR "[INFO] Making '$pi_dir'\n";
  mkdir $pi_dir;
}

## open vcf fofn
open (my $FOFN, $vcf_fofn) or die $!;

## cycle thru fofn
while (my $vcf_file = <$FOFN>) {
  chomp $vcf_file;
  print STDERR "[INFO] > File: '$vcf_file'\n";

  ## iterate thru pops
  foreach my $pop (nsort keys %pops) {
    print STDERR "[INFO] -> Population: '$pop'\n";
    my $sum_pi = 0;

    ## get stats on SNPs in gene region, per population
    open (my $STATS, "bcftools view -R $regions_dir/$gene.regions.txt -S $samples_path/$pop.txt --min-ac=1 --no-update $vcf_file | bcftools stats - | grep \"\^SN\" |") or die $!;
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
      if ( system("bcftools view -R $regions_dir/$gene.regions.txt -S $samples_path/$pop.txt --min-ac=1 --no-update $vcf_file | vcftools --vcf - --out $gene.$pop --site-pi --stdout >$pi_dir/$gene.$pop.sites.pi 2>/dev/null") != 0 ) {
        print STDERR "[INFO] --> Problem with vcftools command!\n";
      } else {
        print STDERR "[INFO] --> Ran vcftools successfully\n";
        open (my $SITESPI, "cut -f3 $pi_dir/$gene.$pop.sites.pi |") or die $!;
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
      $RESULTS{$gene}{$pop}{pi} = 0;
    }
  }
}
close $GENES;

# print Dumper(\%RESULTS);

## print results to file
open (my $RESULTS, ">".$outprefix."_RESULTS.tab") or die $!;
print $RESULTS join ("\t",
  "GENE",
  "LENGTH",
  "POP",
  "N_SAMPLES",
  "N_SNPS",
  "N_MULTIALLELIC_SNPS",
  "N_MNPS",
  "N_INDELS",
  "NUC_DIVERSITY")
  ."\n";

foreach my $k1 (nsort keys %RESULTS) {
  foreach my $k2 (nsort keys %{ $RESULTS{$k1} }) {
    print $RESULTS join ("\t", $k1, $gene_lengths{$k1}, $k2,
      $RESULTS{$k1}{$k2}{num_samples},
      $RESULTS{$k1}{$k2}{num_snps},
      $RESULTS{$k1}{$k2}{num_multiallelic},
      $RESULTS{$k1}{$k2}{num_mnps},
      $RESULTS{$k1}{$k2}{num_indels},
      $RESULTS{$k1}{$k2}{pi} )
      ."\n";

  }
}
close $RESULTS;
print STDERR "[INFO] Finished ".`date`."\n";

__END__
