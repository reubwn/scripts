#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
SYNOPSIS
  Calculates average of a window, e.g. for averaging coverage across a per-site (-d) genomeCov file.

OPTIONS
  -i|--in        [FILE]: infile
  -o|--out       [FILE]: outfile (default = <infile>.<windowsize>.txt)
  -w|--window    [INT] : window size to average across (default=1000)
  -v|--vcf       [FILE]: VCF file input to calculate SNP density across the same window
  -z|--gzip            : input file is gzipped (default=F)
  -f|--flatten         : condense output to one line per window (default=F)
  -p|--pseudochr       : also print 1..n across all scaffolds (default=F)
  -h|--help            : this message

USAGE
  window_cov.pl -i genomeCov.txt -w 1000 -o genomeCov.1000.txt
\n";

my ($infile, $outfile, $vcffile, $gzip, $flatten, $pseudochr, $help);
my $window = 1000;

GetOptions (
  'in|i=s'      => \$infile,
  'out|o:s'     => \$outfile,
  'window|w:i'  => \$window,
  'vcf|v:s'     => \$vcffile,
  'z|gzip'      => \$gzip,
  'flatten|f'   => \$flatten,
  'pseudochr|p' => \$pseudochr,
  'help|h'      => \$help
);

die $usage if $help;
die $usage unless ($infile && $window);
print STDERR "[INFO] Window size: $window\n";
print STDERR "[INFO] Flatten set to TRUE\n" if $flatten;
print STDERR "[INFO] Input is gzipped\n" if $gzip;
unless ($outfile) { $outfile = "$infile.$window.txt" };
my ($IN, $sum);
my (%scaff_lengths, %seen, %v);
my (@string, @averages);
if ($gzip) {
  open ($IN, "zcat $infile |") or die "Cannot open $infile: $!\n";
  $outfile =~ s/\.gz//;
} else {
  open ($IN, $infile) or die "Cannot open $infile: $!\n";
}
open (my $OUT, ">$outfile") or die "Cannot open $outfile: $!\n";
if ($vcffile) {
  print STDERR "[INFO] Analysing VCF file: $vcffile\n";
  my (%e, $cwin);
  open (my $VCF, $vcffile) or die "Cannot open file $vcffile: $!\n";
  while (<$VCF>) {
    chomp;
    my @F = split(/\s+/, $_);
    unless ((/^#/) or (length($F[3])>1)) { ##skip comment lines and any variant that is len>1 (ie not a SNP)
      if (exists($e{$F[0]})){ ##has scaffold been seen before
        if ($F[1] <= ($cwin+$window)) { ##if SNP is in current window range
          $v{$F[0]}{($cwin+$window)}++; ## increment N observed SNPs in that window for that scaffold as HoH
        } else {
          $cwin += $window; ##otherwise add $window to curr window
          if ($F[1] <= ($cwin+$window)) { ##then analyse first SNP of new window
            $v{$F[0]}{($cwin+$window)}++;
          } else {
            $v{$F[0]}{($cwin+$window)} = 0; ##windows with no SNPs get 0
          }
        }
      } else { ##if next scaffold
        $cwin = 0; ##current window is reset
        if ($F[1] <= ($cwin+$window)) { ##and test again
          $v{$F[0]}{($cwin+$window)}++;
        } else {
          $v{$F[0]}{($cwin+$window)} = 0;
        }
        $e{$F[0]}=(); ##seen
      }
    }
  }
  print $OUT "CHROM\tWINDOW\tCOV\tNSNPS\tDENSITY\n";
  close $VCF;
} else {
  print $OUT "CHROM\tWINDOW\tCOV\n";
}

## process genomeCov file
while (<$IN>) {
  chomp;
  push (@string, $_); ##for later printing
  my @F = split(/\s+/, $_);
  $scaff_lengths{$F[0]}++; ##get length of each scaffold
  print STDERR "\r[INFO] Analysing scaffold $F[0]..." if ($. % 1000 == 0); $|=1;

  if ($scaff_lengths{$F[0]} % $window == 0) { ##push if window is reached
    $sum += $F[2]; ##add last coverage
    for (1..$window) {
      push (@averages, ($sum/$window)) ##push to array of windowsize, so each base gets the average value
    }
    if ($flatten) {
      if ($vcffile) {
        print $OUT join("\t", $F[0], $F[1], ($sum/$window), $v{$F[0]}{$F[1]}, ($v{$F[0]}{$F[1]}/$window), "\n")
      } else {
        print $OUT join("\t", $F[0], $F[1], ($sum/$window),"\n");
      }
    }
    $sum = 0; ##reset
   #--0;'''''''''''''''''v  b v hb jb nb   vv                   v cfv hfgyh                                        c  µ ── ·─·─ n baby mamie's line of code :-)
  } else {
    if (eof) { ## special case for eof
      $sum += $F[2];
      for (1..($scaff_lengths{$F[0]} % $window)) {
        push (@averages, ($sum/($scaff_lengths{$F[0]} % $window)));
      }
      if ($flatten) {
        if ($vcffile) {
          print $OUT join("\t", $F[0], $F[1], ($sum/($scaff_lengths{$F[0]} % $window)), $v{$F[0]}{$F[1]}, ($v{$F[0]}{$F[1]}/($scaff_lengths{$F[0]} % $window)), "\n")
        } else {
          print $OUT join("\t", $F[0], $F[1], ($sum/($scaff_lengths{$F[0]} % $window)),"\n");
        }
      }
    }
    $sum += $F[2];
  }
  $seen{$F[0]}=();
}
close $IN;
print STDERR "\n[INFO] Finished averaging...\n[INFO] Printing to $outfile\n";

##print to reannotated file, each site will get average coverage across the specified window
my $n = 1;
unless ($flatten) {
  for my $j (0..$#string) {
    if ($pseudochr) {
      print $OUT join("\t", $string[$j], $n, $averages[$j], "\n");
      $n++;
    } else {
      print $OUT join("\t", $string[$j], $averages[$j], "\n");
    }
  }
}
close $OUT;
print STDERR "[INFO] Finished\n\n";

__END__
