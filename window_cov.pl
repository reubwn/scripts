#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
SYNOPSIS
  Calculates average of a window, e.g. for averaging coverage across a per-site (-d) genomeCov file.
  Generate a per-site coverage file using bedtools:
    >> genomeCoverageBed -ibam <BAMFILE> -g <GENOMEFILE> -d > <BAMFILE>.genomeCov.d

  Can also count the number of SNPs falling within each window if a VCF is provided via the -v flag;
  scaffold names in genomeCov file and VCF file must match up, obviously.

OPTIONS
  -i|--in        [FILE]: infile
  -o|--out       [FILE]: outfile (default = <infile>.<windowsize>.txt)
  -w|--window    [INT] : window size to average across (default=1000)
  -m|--maxcov    [INT] : skip windows with a coverage > this value (default=5000)
  -v|--vcf       [FILE]: VCF file input to calculate SNP density across the same window
  -z|--gzip            : input file is gzipped (default=F)
  -k|--skip            : skip windows < windowsize (ie tailends/small scaffolds) (default=F)
  -h|--help            : this message

USAGE
  window_cov.pl -i genomeCov.txt -w 1000 -o genomeCov.1000.txt
\n";

my ($infile, $outfile, $vcffile, $gzip, $skip, $help);
my $window = 1000;
my $maxcov = 5000;

GetOptions (
  'i|in=s'      => \$infile,
  'o|out:s'     => \$outfile,
  'w|window:i'  => \$window,
  'm|maxcov:i'  => \$maxcov,
  'v|vcf:s'     => \$vcffile,
  'z|gzip'      => \$gzip,
  'k|skip'   => \$skip,
  'h|help'      => \$help
);

die $usage if $help;
die $usage unless ($infile && $window);
print STDERR "[INFO] Window size: $window\n";
print STDERR "[INFO] Input is gzipped\n" if $gzip;
unless ($outfile) { $outfile = "$infile.$window.txt" };

my ($cumulative_coverage_sum, $prev_cumulative_coverage_sum);
my (%cumulative_scaffold_length, %seen, %v);
my (@E);

my $IN;
if ($gzip) {
  open ($IN, "zcat $infile |") or die "[ERROR] Cannot open $infile: $!\n";
  $outfile =~ s/\.gz//;
} else {
  open ($IN, $infile) or die "[ERROR] Cannot open $infile: $!\n";
}
open (my $OUT, ">$outfile") or die "[ERROR] Cannot open $outfile: $!\n";
if ($vcffile) {
  print STDERR "[INFO] Analysing VCF file: $vcffile\n";
  my (%e, $cwin);
  open (my $VCF, $vcffile) or die "[ERROR] Cannot open file $vcffile: $!\n";
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
LINE: while (<$IN>) {
  chomp;
  my @F = split(/\s+/, $_);
  print STDERR "\r[INFO] Analysing scaffold $F[0]..." if ($. % 5000 == 0); $|=1; ##just print info every 5000 bases

  ## add coverage at current site to cumulative_coverage_sum:
  $cumulative_coverage_sum += $F[2];

  ## print if window size is reached:
  if ($. == 1) {
    ## empty condition to capture the first line
  } elsif ( $F[1] % $window == 0) {
    if ($vcffile) {
      my $nsnps;
      if (exists($v{$F[0]}{$F[1]})) { $nsnps = $v{$F[0]}{$F[1]} } else { $nsnps = 0 };
      print $OUT join("\t", $F[0], $F[1], ($cumulative_coverage_sum/$window), $nsnps, ($nsnps/$window), "\n");
    } else {
      print $OUT join("\t", $F[0], $F[1], ($cumulative_coverage_sum/$window),"\n");
    }

    ## reset cumsum to 0:
    $cumulative_coverage_sum = 0;

  ## OR its a new scaffold (not seen already) OR eof:
  } elsif ( (!exists($seen{$F[0]})) || (eof) ) {
    unless ($E[1] == $window) { ## avoids double processing and error for cases where scaffold length is exactly == window size
      if ($vcffile) {
        my $nsnps;
        if (exists($v{$E[0]}{$E[1]})) { $nsnps = $v{$E[0]}{$E[1]} } else { $nsnps = 0 }; ##WARNING: THIS MIGHT NEVER WORK.....!!!!!
        print $OUT join("\t", $E[0], $E[1], ($prev_cumulative_coverage_sum/($E[1] % $window)), $nsnps, ($nsnps/($E[1] % $window)), "\n"); ##WARNING: NSNPS and DENSITY will always be zero for this last window.... IS THIS OK
      } else {
        print $OUT join("\t", $E[0], $E[1], ($prev_cumulative_coverage_sum/($E[1] % $window)),"\n");
      }

      ## reset cumsum to 0 then begin again:
      $cumulative_coverage_sum = 0;
      $cumulative_coverage_sum += $F[2];
    }
  }

  @E = @F; ##store current line in variable accessible to next iteration
  $prev_cumulative_coverage_sum = $cumulative_coverage_sum; ##store current cumsum in variable accessible to next iteration
  $seen{$F[0]}=(); ##see the scaffold

  # if ($cumulative_scaffold_length{$F[0]} % $window == 0) { ##push or print if window is reached
  #   $sum += $F[2]; ##add last coverage
  #   for (1..$window) {
  #     push (@averages, ($sum/$window)) ##push to array of windowsize, so each base gets the average value
  #   }
  #   if ($flatten) {
  #     if ($vcffile) {
  #       my $nsnps;
  #       if (exists($v{$F[0]}{$F[1]})) { $nsnps = $v{$F[0]}{$F[1]} } else { $nsnps = 0 };
  #       print $OUT join("\t", $F[0], $F[1], ($sum/$window), $nsnps, ($nsnps/$window), "\n")
  #     } else {
  #       print $OUT join("\t", $F[0], $F[1], ($sum/$window),"\n");
  #     }
  #   }
  #   $sum = 0; ##reset
  #  #--0;'''''''''''''''''v  b v hb jb nb   vv                   v cfv hfgyh                                        c  µ ── ·─·─ n baby mamie's line of code :-)
  # } else {
  #   if (eof) { ## special case for eof
  #     $sum += $F[2];
  #     for (1..($cumulative_scaffold_length{$F[0]} % $window)) {
  #       push (@averages, ($sum/($cumulative_scaffold_length{$F[0]} % $window)));
  #     }
  #     if ($flatten) {
  #       if ($vcffile) {
  #         my $nsnps;
  #         if (exists($v{$F[0]}{$F[1]})) { $nsnps = $v{$F[0]}{$F[1]} } else { $nsnps = 0 };
  #         print $OUT join("\t", $F[0], $F[1], ($sum/($cumulative_scaffold_length{$F[0]} % $window)), $nsnps, ($nsnps/($cumulative_scaffold_length{$F[0]} % $window)), "\n")
  #       } else {
  #         print $OUT join("\t", $F[0], $F[1], ($sum/($cumulative_scaffold_length{$F[0]} % $window)),"\n");
  #       }
  #     }
  #   }
  #   $sum += $F[2];
  # }
  # $seen{$F[0]}=();
}
close $IN;
print STDERR "\n[INFO] Finished averaging...\n[INFO] Printing to $outfile\n";

# ##print to reannotated file, each site will get average coverage across the specified window
# my $n = 1;
# unless ($flatten) {
#   for my $j (0..$#string) {
#     if ($pseudochr) {
#       print $OUT join("\t", $string[$j], $n, $averages[$j], "\n");
#       $n++;
#     } else {
#       print $OUT join("\t", $string[$j], $averages[$j], "\n");
#     }
#   }
# }
close $OUT;
print STDERR "[INFO] Finished on ".`date`."\n\n";

__END__
