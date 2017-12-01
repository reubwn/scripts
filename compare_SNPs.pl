#!/usr/bin/env perl

## author: reubwn Nov 2017

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;

my $usage = "
OPTIONS:
  -1|--file1 [FILE]
  -2|--file2 [FILE]
  -h|--help
\n";

my ($file1,$file2,$help);
my $outfile = "compare";
my $processed = 0;

GetOptions (
  '1|file1=s' => \$file1,
  '2|file2=s' => \$file2,
  'o|out:s'   => \$outfile,
  'h|help'    => \$help,
);

die $usage if $help;
die $usage unless ($file1 && $file2);

my (%h1,%h2,%intersect,%uniq1,%uniq2);
my ($mismatch,$uniq1,$uniq2) = (0,0,0);

## parse SNPs in file1/2:
open (my $FILE1, $file1) or die $!;
while (<$FILE1>) {
  chomp;
  my @F = split /\s+/;
  $h1{"$F[0].$F[1]"} = { 'chrom' => $F[0], 'pos' => $F[1], 'ref' => $F[2], 'alt' => $F[3], 'TC' => $F[5], 'TR' => $F[6], 'MAF' => $F[8] }; ##key= pos; val= %{chrom...}
}
close $FILE1;
open (my $FILE2, $file2) or die $!;
while (<$FILE2>) {
  chomp;
  my @F = split /\s+/;
  $h2{"$F[0].$F[1]"} = { 'chrom' => $F[0], 'pos' => $F[1], 'ref' => $F[2], 'alt' => $F[3], 'TC' => $F[5], 'TR' => $F[6], 'MAF' => $F[8] }; ##key= pos; val= %{chrom...}
}
close $FILE2;

open (my $INTERSECT, ">$outfile.intersect") or die $!;
open (my $UNIQ1, ">$outfile.uniq.1") or die $!;
open (my $UNIQ2, ">$outfile.uniq.2") or die $!;
foreach (nsort keys %h1) {
  ## SNP exists in same position on same CHROM:
  if ( (exists($h2{$_})) and ($h1{$_}{chrom} eq $h2{$_}{chrom}) ) {
    ## check REF and ALT alleles are also the same:
    unless ( ($h1{$_}{ref} eq $h2{$_}{ref}) and ($h1{$_}{alt} eq $h2{$_}{alt}) ) {
      print STDERR "[INFO] SNP at position $h1{$_}{chrom}:$h1{$_}{pos} has $h1{$_}{ref}/$h1{$_}{alt} in $file1 but $h2{$_}{ref}/$h2{$_}{alt} in $file2\n";
      $mismatch++;
    } else {
      print $INTERSECT join (
        "\t",
        $h1{$_}{chrom},
        $h1{$_}{pos},
        $h1{$_}{ref},
        $h1{$_}{alt},
        $h1{$_}{TC},
        $h1{$_}{TR},
        $h1{$_}{MAF},
        $h2{$_}{chrom},
        $h2{$_}{pos},
        $h2{$_}{ref},
        $h2{$_}{alt},
        $h2{$_}{TC},
        $h2{$_}{TR},
        $h2{$_}{MAF},
        "\n"
      );
      $intersect{$_}++;
    }
  } else {
    print $UNIQ1 join (
      "\t",
      $h1{$_}{chrom},
      $h1{$_}{pos},
      $h1{$_}{ref},
      $h1{$_}{alt},
      $h1{$_}{TC},
      $h1{$_}{TR},
      $h1{$_}{MAF},
      "\n"
    );
  }
}
close $INTERSECT;
close $UNIQ1;

foreach (keys %h2) {
  unless (exists($intersect{$_})) {
    print $UNIQ2 join (
      "\t",
      $h2{$_}{chrom},
      $h2{$_}{pos},
      $h2{$_}{ref},
      $h2{$_}{alt},
      $h2{$_}{TC},
      $h2{$_}{TR},
      $h2{$_}{MAF},
      "\n"
    );
  }
}
close $UNIQ2;

print STDERR "[INFO] # SNPs in $file1: ".commify(scalar(keys %h1))."\n";
print STDERR "[INFO] # SNPs in $file2: ".commify(scalar(keys %h2))."\n";
print STDERR "[INFO] # SNPs common to both files: ".commify(scalar(keys %intersect))."\n";
print STDERR "[INFO] # SNPs at same position but mismatched states: ".commify($mismatch)."\n";
print STDERR "[INFO]   % SNPs $file1: ".percentage(scalar(keys %intersect),scalar(keys %h1))."\n";
print STDERR "[INFO]   % SNPs $file2: ".percentage(scalar(keys %intersect),scalar(keys %h2))."\n";

################### SUBS

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
