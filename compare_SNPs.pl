#!/usr/bin/env perl

## author: reubwn Nov 2017

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
OPTIONS:
  -1|--file1 [FILE] : SNP freqs file 1
  -2|--file2 [FILE] : SNP freqs file 2
  -f|--fasta [FILE] : reference fasta sequence file (optional)
  -o|--out   [STR]  : output prefix (default 'compare_SNP')
  -h|--help
\n";

my ($file1,$file2,$fasta,$help);
my $outfile = "compare_SNP";
my $processed = 0;

GetOptions (
  '1|file1=s' => \$file1,
  '2|file2=s' => \$file2,
  'f|fasta=s' => \$fasta,
  'o|out:s'   => \$outfile,
  'h|help'    => \$help,
);

die $usage if $help;
die $usage unless ($file1 && $file2);

my (%seqlengths,%h1,%h2,%intersect,%uniq1,%uniq2);
my ($mismatch,$uniq1,$uniq2) = (0,0,0);


print STDERR "[INFO] Getting sequence lengths from: $fasta...\n";
my $in = Bio::SeqIO->new( -file => $fasta, -format => 'fasta' );
while ( my $seqobj = $in -> next_seq() ) {
  $seqlengths{$seqobj->display_id()} = $seqobj->length();
}

## parse SNPs in file1/2:
open (my $FILE1, $file1) or die $!;
while (<$FILE1>) {
  chomp;
  my @F = split /\s+/;
  $h1{"$F[0].$F[1]"} = { 'chrom' => $F[0], 'pos' => $F[1], 'ref' => $F[2], 'alt' => $F[3], 'TC' => $F[5], 'TR' => $F[6], 'MAF' => $F[9] }; ##key= pos; val= %{chrom...}
}
close $FILE1;
open (my $FILE2, $file2) or die $!;
while (<$FILE2>) {
  chomp;
  my @F = split /\s+/;
  $h2{"$F[0].$F[1]"} = { 'chrom' => $F[0], 'pos' => $F[1], 'ref' => $F[2], 'alt' => $F[3], 'TC' => $F[5], 'TR' => $F[6], 'MAF' => $F[9] }; ##key= pos; val= %{chrom...}
}
close $FILE2;

open (my $INTERSECT, ">$outfile.intersect") or die $!;
open (my $MISMATCH, ">$outfile.mismatch") or die $!;
open (my $UNIQ1, ">$outfile.uniq.1") or die $!;
open (my $UNIQ2, ">$outfile.uniq.2") or die $!;
open (my $INFO, ">$outfile.summary") or die $!;
foreach (sort {ncmp($h1{$a}{chrom},$h1{$b}{chrom})} keys %h1) {
  ## SNP exists in same position on same CHROM:
  if ( (exists($h2{$_})) and ($h1{$_}{chrom} eq $h2{$_}{chrom}) ) {
    ## check REF and ALT alleles are also the same:
    unless ( ($h1{$_}{ref} eq $h2{$_}{ref}) and ($h1{$_}{alt} eq $h2{$_}{alt}) ) {
      print $MISMATCH join (
        "\t",
        $h1{$_}{chrom},
        $seqlengths{$h1{$_}{chrom}},
        $h1{$_}{pos},
        $h1{$_}{ref},
        $h1{$_}{alt},
        $h1{$_}{TC},
        $h1{$_}{TR},
        $h1{$_}{MAF},
        $h2{$_}{chrom},
        $seqlengths{$h2{$_}{chrom}},
        $h2{$_}{pos},
        $h2{$_}{ref},
        $h2{$_}{alt},
        $h2{$_}{TC},
        $h2{$_}{TR},
        $h2{$_}{MAF},
        "\n"
      );
      $mismatch++;
    } else {
      print $INTERSECT join (
        "\t",
        $h1{$_}{chrom},
        $seqlengths{$h1{$_}{chrom}},
        $h1{$_}{pos},
        $h1{$_}{ref},
        $h1{$_}{alt},
        $h1{$_}{TC},
        $h1{$_}{TR},
        $h1{$_}{MAF},
        $h2{$_}{chrom},
        $seqlengths{$h2{$_}{chrom}},
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
      $seqlengths{$h1{$_}{chrom}},
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
close $MISMATCH;
close $UNIQ1;

foreach (sort {ncmp($h2{$a}{chrom},$h2{$b}{chrom})} keys %h2) {
  unless (exists($intersect{$_})) {
    print $UNIQ2 join (
      "\t",
      $h2{$_}{chrom},
      $seqlengths{$h2{$_}{chrom}},
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
print STDERR "[INFO]   % SNPs $file1: ".percentage(scalar(keys %intersect),scalar(keys %h1))."\n";
print STDERR "[INFO]   % SNPs $file2: ".percentage(scalar(keys %intersect),scalar(keys %h2))."\n";
print STDERR "[INFO] # SNPs at same position but mismatched states: ".commify($mismatch)."\n";
print STDERR "[INFO] Finished on ".`date`."\n";

print $INFO "[INFO] # SNPs in $file1: ".commify(scalar(keys %h1))."\n";
print $INFO "[INFO] # SNPs in $file2: ".commify(scalar(keys %h2))."\n";
print $INFO "[INFO] # SNPs common to both files: ".commify(scalar(keys %intersect))."\n";
print $INFO "[INFO]   % SNPs $file1: ".percentage(scalar(keys %intersect),scalar(keys %h1))."\n";
print $INFO "[INFO]   % SNPs $file2: ".percentage(scalar(keys %intersect),scalar(keys %h2))."\n";
print $INFO "[INFO] # SNPs at same position but mismatched states: ".commify($mismatch)."\n";
print $INFO "[INFO] Finished on ".`date`."\n";
close $INFO;

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
