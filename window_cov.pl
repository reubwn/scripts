#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
SYNOPSIS
  Calculates average of a window, e.g. for averaging coverage across a per-site (-d) genomeCov file.

USAGE
  window_cov.pl -i genomeCov.txt -w 1000 -o genomeCov.1kb.txt
\n";

my ($infile, $help);
my $window = 1000;
my $outfile = "$infile.$window.txt";

GetOptions (
  'in|i=s'     => \$infile,
  'out|o:s'    => \$outfile,
  'window|w:i' => \$window,
  'help|h'     => \$help
);

die $usage if $help;
die $usage unless ($infile && $window);
print STDERR "[INFO] Window size: $window\n";

my $sum;
my (%scaff_lengths, %seen);
my (@string,@averages);
open (my $IN, $infile) or die "Cannot open $infile: $!\n";

while (<$IN>) {
  chomp;
  push (@string, $_); ##for later printing
  my @F = split(/\s+/, $_);
  $scaff_lengths{$F[0]}++; ##get length of each scaffold
  print STDERR "\r[INFO] Scaffold $F[0]...";

  if ($scaff_lengths{$F[0]} % $window == 0) { ## push if window is reached
    $sum += $F[2]; ##add last coverage
    for (1..$window) {
      push (@averages, ($sum/$window)) ## push to array of windowsize, so each base gets the average value
    }
    $sum = 0; ##reset

  } else {
    if (!exists($seen{$F[0]})) { ##special case when scaffold name changes
      unless ($.==1) { ## but not on first instance
        for (1..($scaff_lengths{"RMAC00001"} % $window)) {
          push (@averages, ($sum/($scaff_lengths{"RMAC00001"} % $window)));
        }
        $sum = 0;
      }
    }

    if (eof) { ## special case for eof
      $sum += $F[2];
      for (1..($scaff_lengths{$F[0]} % $window)) {
        push (@averages, ($sum/($scaff_lengths{$F[0]} % $window)));
      }
    }

    $sum += $F[2];
  }
  $seen{$F[0]}=();
}
close $IN;
print STDERR "\n[INFO] Finished averaging...\n[INFO] Printing to $outfile\n";

##print to reannotated file, each site will get average coverage across the specified window
open (my $OUT, ">$outfile") or die "Cannot open $outfile: $!\n";
for my $j (0..$#string) {
  print $OUT join("\t", $string[$j], $averages[$j], "\n");
}
close $OUT;
print STDERR "[INFO] Finished\n\n";

__END__
