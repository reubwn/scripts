#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw /sum/;

my $usage = "
SYNOPSIS:
  Calculates 'collinearity' score based on the number of collinear genes divided by the total number of genes within that defined block.
  Takes the collinearity file and the 'gff' file used in MCScanX analyses.

OUTPUT:
  Prints to a file 'xyz.collinearity.score'; prints score for each block plus an average.

OPTIONS:
  -c|--collinearity [FILE] : collinearity file from MCScanX
  -g|--gff          [FILE] : modified gff file from MCScanX
  -k|--kaks                : parse collinearity file to get average ka & ks per block
  -h|--help                : print this message

USAGE:
  >> calculate_collinarity_metric.pl -c xyz.collinearity -g xyz.gff
\n";

my ($collinearity, $gff, $kaks, $help);

GetOptions (
  'collinearity|c=s' => \$collinearity,
  'gff|g=s'          => \$gff,
  'kaks|k'           => \$kaks,
  'help|h'           => \$help,
);

die $usage if $help;
die $usage unless ($collinearity && $gff);

my (%blocks,$orientation);
open (my $COL, $collinearity) or die $!;
while (<$COL>) {
  chomp;
  my ($aln_number);
  if ($_ =~ m/^#/) {
    if ($_ =~ m/(plus|minus)$/) { ## get strand orientation of block 2
      $orientation = $1;
      next;
    } else {
      next;
    }
  }
  $_ =~ s/^\s+|\s+$//g; ##remove leading and trailing whitespaces

  my @F = split (m/\s+/, $_);
  if ($F[0]=~m/^\d+\-\d+\:$/) { ## sometimes columns not formatted properly... :/
    my @a = split (m/\-/, $F[0]);
    push @{ $blocks{$a[0]}{'block1'} }, $F[1];
    push @{ $blocks{$a[0]}{'block2'} }, $F[2];
    $aln_number = $a[0];
  } else {
    $F[0] =~ s/\-//;
    push @{ $blocks{$F[0]}{'block1'} }, $F[2];
    push @{ $blocks{$F[0]}{'block2'} }, $F[3];
    $aln_number = $F[0];
    #$aln_number =~ s/\-//;
  }

  ## dump genes and plus/minus info into %blocks
  #push @{ $blocks{$aln_number}{'block1'} }, $F[2];
  #push @{ $blocks{$aln_number}{'block2'} }, $F[3];
  $blocks{$aln_number}{ 'orientation'} = $orientation;

  if ($kaks) {
    push @{ $blocks{$aln_number}{'ks'} }, $F[-1]; ## ks is in final column
    push @{ $blocks{$aln_number}{'ka'} }, $F[-2]; ## ka is in second to last column
  }
}
close $COL;

open (my $OUT, ">$collinearity.score") or die $!;
if ($kaks) {
  print $OUT join "\t", "block_num", "collinear_genes", "total_genes1", "total_genes2", "orientation", "score_block1", "score_block2", "score_avg", "ka_avg", "ks_avg", "\n";
} else {
  print $OUT join "\t", "block_num", "collinear_genes", "total_genes1", "total_genes2", "orientation", "score_block1", "score_block2", "score_avg", "\n";
}

foreach (sort {$a<=>$b} keys %blocks) {

  ## get orientation of block2
  my $orientation = $blocks{$_}{'orientation'};

  ## get genes of block1
  my @block1_genes = @{ $blocks{$_}{'block1'} };
  my $bl1_start = shift @block1_genes;
  my $bl1_end   = pop @block1_genes;
  my $bl1_length = `perl -e 'while (<>){print if (/\t\Q$bl1_start\E\t/../\t\Q$bl1_end\E\t/);}' $gff | wc -l`;
  chomp ($bl1_length);
  my $score_block1 = sprintf("%.5f",(scalar(@block1_genes)/$bl1_length));

  ## get genes of block2
  my @block2_genes = @{ $blocks{$_}{'block2'} };
  my $bl2_start = shift @block2_genes;
  my $bl2_end   = pop @block2_genes;
  my $bl2_length;
  ## if block 2 is in minus orientation, need to reverse the search!
  if ($orientation eq "plus") {
    $bl2_length = `perl -e 'while (<>){print if (/\t\Q$bl2_start\E\t/../\t\Q$bl2_end\E\t/);}' $gff | wc -l`;
  } elsif ($orientation eq "minus") {
    $bl2_length = `perl -e 'while (<>){print if (/\t\Q$bl2_end\E\t/../\t\Q$bl2_start\E\t/);}' $gff | wc -l`;
  } else {
    die "\nUnknown strand orientation for block 2: $orientation\n\n";
  }
  chomp ($bl2_length);
  print STDOUT "\rCalculating scores for block: $_";
  my $score_block2 = sprintf("%.5f",(scalar(@block2_genes)/$bl2_length));

  ## get kaks values if present
  if ($kaks) {
    my @ka = grep {$_ >= 0} @{ $blocks{$_}{'ka'} }; ## exclude negative values from calculation;
    my @ks = grep {$_ >= 0} @{ $blocks{$_}{'ks'} }; ## these are "-2", output when ka/ks cannot be calulated for some reason

    print $OUT join "\t",
      $_,
      scalar(@block1_genes),
      $bl1_length,
      $bl2_length,
      $orientation,
      $score_block1,
      $score_block2,
      sprintf("%.5f",(($score_block1+$score_block2)/2)),
      sprintf("%.5f",(avg(@ka))),
      sprintf("%.5f",(avg(@ks))),
      "\n";
  } else {
    print $OUT join "\t",
      $_,
      scalar(@block1_genes),
      $bl1_length,
      $bl2_length,
      $orientation,
      $score_block1,
      $score_block2,
      sprintf("%.5f",(($score_block1+$score_block2)/2)),
      "\n";
  }
}
close $OUT;

sub avg {
  return sum(@_)/@_;
}

__END__
