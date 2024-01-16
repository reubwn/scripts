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
  Sample filtering options work on gene-by-gene basis.

GENERAL OPTIONS [*required]
  -v|--vcf      *[FILE]  : VCF/BCF file
  -g|--genes    *[FILE]  : TXT list of gene IDs to be analysed (linking to the 'mRNA' field of the GFF)
  -G|--gff      *[FILE]  : GFF annotation file
  -p|--pops     *[FILE]  : TXT list of population/site IDs to subset VCF
  -s|--samples  *[PATH]  : path/to/dir containing samples files that link sample ID to population grouping
  -o|--out       [STR]   : outfiles prefix ('pi')
  -h|--help              : print this message

SAMPLE FILTERING OPTIONS
  -e|--het       [FLOAT] : proportion of heterozygous calls per sample allowed (default = 0.5; set to 1 to disable)
  -m|--missing   [FLOAT] : proportion of missing data per sample (default = 0.5; set to 1 to disable)
  -n|--samples_N [INT]   : minimum number of samples per population after applying above filters (default = 5)

\n";

my ($vcf_file, $genes_file, $gff_file, $pops_file, $samples_path, $help);
my $outprefix = "pi";
my $het_threshold = 0.5;
my $missing_threshold = 0.5;
my $samples_N_threshold = 5;

GetOptions (
  'v|vcf=s'     => \$vcf_file,
  'g|genes=s'   => \$genes_file,
  'G|gff=s'     => \$gff_file,
  'p|pops=s'    => \$pops_file,
  's|samples=s' => \$samples_path,
  'e|het:f'     => \$het_threshold,
  'm|missing:f' => \$missing_threshold,
  'n|samples_N:i' => \$samples_N_threshold,
  'o|out:s'     => \$outprefix,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $vcf_file && $genes_file && $pops_file && $samples_path );

## RESULTS hash
my %RESULTS;
my %gene_lengths;

## parse populations files
open (my $POPS, $pops_file) or die "Problem with '$pops_file': $!";
my %pops;
while (my $pop = <$POPS>) {
  chomp $pop;
  $pops{$pop}++;
}
close $POPS;
print STDERR "[INFO] VCF file: '$vcf_file'\n";
print STDERR "[INFO] GFF file: '$gff_file'\n";
print STDERR "[INFO] Genes file: '$genes_file'\n";
print STDERR "[INFO] Populations file: '$pops_file'\n";
print STDERR "[INFO] Number of populations in '$pops_file': ".scalar(keys %pops)."\n";
print STDERR "[INFO] Threshold for missing data = $missing_threshold\n";
print STDERR "[INFO] Threshold for heterozygous calls = $het_threshold\n";
print STDERR "[INFO] Threshold samples size (after filtering) = $samples_N_threshold\n\n";

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

## open genes file
open (my $GENES, $genes_file) or die "Problem with '$genes_file': $!";

## grep gene IDs from GFF to generate regions.txt file
while (my $gene = <$GENES>) {
  chomp $gene;
  print STDERR "\n[####]\n[INFO] Gene: '$gene'\n[####]\n";
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
  POP: foreach my $pop (nsort keys %pops) {
    my $sum_pi = 0;
    ## samples to exclude based on missing data and/or proportion of het calls
    my @excluded_samples_het;
    my @excluded_samples_missing;

    ## get stats on SNPs in gene region, per population
    open (my $STATS, "bcftools view -R $regions_dir/$gene.regions.txt -S $samples_path/$pop.txt -Ou $vcf_file | bcftools view -a -c1 -Ou | bcftools view -i 'TYPE=\"snp\"' -Ou | bcftools stats -s- |") or die $!;
    while (my $line = <$STATS>) {
      chomp $line;
      if ($line =~ m/^SN/) { ## summary numbers block
        my @F = split (/\s+/, $line);
        $RESULTS{$gene}{$pop}{num_samples} = $F[-1] if $line =~ "number of samples";
        $RESULTS{$gene}{$pop}{num_snps} = $F[-1] if $line =~ "number of SNPs";
        $RESULTS{$gene}{$pop}{num_mnps} = $F[-1] if $line =~ "number of MNPs";
        $RESULTS{$gene}{$pop}{num_indels} = $F[-1] if $line =~ "number of indels";
        $RESULTS{$gene}{$pop}{num_multiallelic} = $F[-1] if $line =~ "number of multiallelic SNP sites";

      } elsif ($line =~ m/^PSC/) { ## per-sample counts block
        if ($RESULTS{$gene}{$pop}{num_snps} > 0) { ## only do for genes with SNPs
          my @F = split (/\s+/, $line);
          ## sum of 4th, 5th, 6th and 14th cols = total num SNPs in subset VCF
          # PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.
          # PSC	[2]id	[3]sample	[4]nRefHom	[5]nNonRefHom	[6]nHets	[7]nTransitions	[8]nTransversions	[9]nIndels	[10]average depth	[11]nSingletons	[12]nHapRef	[13]nHapAlt	[14]nMissing
          push (@excluded_samples_het, $F[2]) if (($F[5]/($F[3]+$F[4]+$F[5]+$F[13])) > $het_threshold);
          push (@excluded_samples_missing, $F[2]) if (($F[13]/($F[3]+$F[4]+$F[5]+$F[13])) > $missing_threshold);
        }
      }
    }
    close $STATS;
    print STDERR "[INFO] Population: '$pop' (N = $RESULTS{$gene}{$pop}{num_samples})\n";

    ## join lists, add info to RESULTS
    my @excluded_all = (@excluded_samples_het, @excluded_samples_missing);
    $RESULTS{$gene}{$pop}{num_excluded} = scalar(@excluded_all);
    $RESULTS{$gene}{$pop}{excluded_samples_het} = join(",",@excluded_samples_het);
    $RESULTS{$gene}{$pop}{excluded_samples_missing} = join(",",@excluded_samples_missing);
    $RESULTS{$gene}{$pop}{excluded_samples_all} = join(",",@excluded_all);
    $RESULTS{$gene}{$pop}{num_samples_final} = ($RESULTS{$gene}{$pop}{num_samples} - $RESULTS{$gene}{$pop}{num_excluded}); ## num samples after excluding some samples based on missing data and/or het calls

    ## iterate to next pop if N < sample_N_threshold after excluding additional samples
    if ( $RESULTS{$gene}{$pop}{num_samples_final} < $samples_N_threshold ) {
      print STDERR "[INFO] Pop '$pop' has < $samples_N_threshold samples after filtering, skipping\n";
      next POP;
    }

    ## execute slightly different commands depending on whether additional sample filtering is required
    if (scalar(@excluded_all) > 0) {
      print STDERR "[INFO] \tTotal samples excluded: ".scalar(@excluded_all)." (".percentage(scalar(@excluded_all),$RESULTS{$gene}{$pop}{num_samples},1).")\n";
      # print STDERR "[INFO] Total samples excluded: ".scalar(@excluded_all)." (MISSING>$missing_threshold = ".scalar(@excluded_samples_missing)."; HET>$het_threshold = ".scalar(@excluded_samples_het).")\n";
      my $exclude_string = join(",", @excluded_all);

      ## command block if samples are to be excluded
      ## skip if no SNPs
      if ( $RESULTS{$gene}{$pop}{num_snps} > 0 ) {
        print STDERR "[INFO] \tNumber SNPs: $RESULTS{$gene}{$pop}{num_snps}\n";
        ## run vcftools --site-pi
        if ( system("bcftools view -Ou -R $regions_dir/$gene.regions.txt -S $samples_path/$pop.txt $vcf_file | bcftools view -Ou -s ^$exclude_string | bcftools view -Ou -a -c1 | bcftools view -Ov -i 'TYPE=\"snp\"' | vcftools --vcf - --out $gene.$pop --site-pi --stdout >$pi_dir/$gene.$pop.sites.pi 2>/dev/null") != 0 ) {
          print STDERR "[INFO] Problem with vcftools command!\n\n";
          die 1;
        } else {
          # print STDERR "[INFO] Ran vcftools successfully\n";
          open (my $SITESPI, "cut -f3 $pi_dir/$gene.$pop.sites.pi |") or die $!;
          while (my $pi = <$SITESPI>) {
            next if $. == 1;
            chomp $pi;
            $sum_pi += $pi;
          }
          close $SITESPI;
          $RESULTS{$gene}{$pop}{pi} = ($sum_pi/$gene_lengths{$gene});
          print STDERR "[INFO] \tNucleotide diversity = $RESULTS{$gene}{$pop}{pi}\n";
        }
      } else {
        print STDERR "[INFO] \tNo variants found\n";
        $RESULTS{$gene}{$pop}{pi} = 0;
      }

    } else {
      ## command block if no samples are to be excluded
      ## skip if no SNPs
      if ( $RESULTS{$gene}{$pop}{num_snps} > 0 ) {
        print STDERR "[INFO] \tNumber SNPs: $RESULTS{$gene}{$pop}{num_snps}\n";
        ## run vcftools --site-pi
        if ( system("bcftools view -Ou -R $regions_dir/$gene.regions.txt -S $samples_path/$pop.txt $vcf_file | bcftools view -Ou -a -c1 | bcftools view -Ov -i 'TYPE=\"snp\"' | vcftools --vcf - --out $gene.$pop --site-pi --stdout >$pi_dir/$gene.$pop.sites.pi 2>/dev/null") != 0 ) {
          print STDERR "[INFO] Problem with vcftools command!\n\n";
          die 1;
        } else {
          # print STDERR "[INFO] Ran vcftools successfully\n";
          open (my $SITESPI, "cut -f3 $pi_dir/$gene.$pop.sites.pi |") or die $!;
          while (my $pi = <$SITESPI>) {
            next if $. == 1;
            chomp $pi;
            $sum_pi += $pi;
          }
          close $SITESPI;
          $RESULTS{$gene}{$pop}{pi} = ($sum_pi/$gene_lengths{$gene});
          print STDERR "[INFO] \tNucleotide diversity = $RESULTS{$gene}{$pop}{pi}\n";
        }
      } else {
        print STDERR "[INFO] \tNo variants found\n";
        $RESULTS{$gene}{$pop}{pi} = 0;
      }
    }
  }
  # print STDERR "\n";
}
close $GENES;

# print Dumper(\%RESULTS);

## print results to file
open (my $RESULTS, ">".$outprefix."_RESULTS.tab") or die $!;
print $RESULTS join ("\t",
  "GENE",
  "LENGTH",
  "POP",
  "N_SAMPLES_INIT",
  "N_EXCLUDED",
  "EXCLUDED_IDS",
  "N_SAMPLES_FINAL",
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
      $RESULTS{$k1}{$k2}{num_excluded},
      $RESULTS{$k1}{$k2}{excluded_samples_all},
      $RESULTS{$k1}{$k2}{num_samples_final},
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

######################### sub-routines

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

sub percentage {
  my $numerator = $_[0];
  my $denominator = $_[1];
  my $places = "\%.2f"; ## default is two decimal places
  if (exists $_[2]){$places = "\%.".$_[2]."f";};
  my $float = (($numerator / $denominator)*100);
  my $rounded = sprintf("$places",$float);
  return "$rounded\%";
}

__END__
