#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use List::Util qw/sum/;

my $usage = "
  Get intergenic distances from a GFF3 file

  USAGE: get_intergenic.pl [options] <GFFFILE>
  OPTIONS:
    -s|--scaffold : prints mean intergenic dist per scaffold
    -h|--help     : this message\n
";

my ($scaffold,$help);
GetOptions (
  's|scaffold' => \$scaffold,
  'h|help'     => \$help
);

die $usage if $help;
die $usage if @ARGV == 0;

my %h;
open (my $GFF, $ARGV[0]) or die $!;
while (<$GFF>) {
  chomp;
  my @F = split (/\s+/, $_);
  if ($F[2] =~ "gene") {
    push ( @{$h{$F[0]}{STARTS}}, $F[3] );
    push ( @{$h{$F[0]}{ENDS}}, $F[4] );
  } else {
    next;
  }
}
close $GFF;

my %scaffold;
foreach my $chrom (nsort keys %h) {
  my @STARTS = @{$h{$chrom}{STARTS}};
  my @ENDS = @{$h{$chrom}{ENDS}};
  #print "$chrom\n@STARTS\n@ENDS\n";
  for my $i (1 .. $#STARTS) {
    my $distance = $STARTS[$i] - $ENDS[($i-1)];
    if ($scaffold) {
      push ( @{$scaffold{$chrom}}, $distance ) unless $distance < 0;
    } else {
      print "$chrom\t$distance\n" unless $distance < 0;
    }
  }
}

if ($scaffold) {
  foreach my $chrom (nsort keys %scaffold) {
    my @a = @{$scaffold{$chrom}};
    my $mean = sum(@a)/scalar(@a);
    print "$chrom\t$mean\n";
  }
}
