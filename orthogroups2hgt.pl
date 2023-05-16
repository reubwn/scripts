#!/usr/bin/env perl

## reubwn Aug 2018

use strict;
use warnings;

use Getopt::Long;
use Term::ANSIColor;
use Sort::Naturally;

my $usage = "
SYNOPSIS
  Paints the Orthogroups.txt file with annotations from the HGT_results file.

OPTIONS
  -i|--orthogroups [FILE] : Orthogroups.txt file from OrthoFinder
  -a|--annot       [FILE] : *HGT_results file(s)
  -u|--hU          [INT]  : hU threshold for defining HGTc (>30)
  -c|--CHS         [INT]  : CHS threshold for defining HGTc (>90\%)
  -s|--simple             : print simple output
  -o|--outfile     [STR]  : output base filename ('inputfilename')
  -h|--help               : print this message
\n";

my ($orthogroupsfile, $annot, $simple, $outfile, $help);
my $hU_threshold = 30;
my $CHS_threshold = 90;

GetOptions (
  'i|orthogroups=s' => \$orthogroupsfile,
  'a|annot=s'       => \$annot,
  'u|hU:i'      => \$hU_threshold,
  'c|CHS:i'     => \$CHS_threshold,
  's|simple'    => \$simple,
  'o|outfile'   => \$outfile,
  'h|help'          => \$help
);

die $usage if $help;
die $usage unless ($orthogroupsfile && $annot);

## outfiles
$outfile = $orthogroupsfile.".annot" unless ($outfile);
open (my $OUT, ">$outfile") or die $!;
open (my $PROP, ">$outfile.proportion") or die $!;

## parse $annot if present
my %annot_hash;
print STDERR "[INFO] Collecting annotations from " . colored($annot, 'white on_blue') . "\n";
open (my $ANNOT, $annot) or die $!;
LINE: while (my $line = <$ANNOT>) {
  next LINE if $line =~ m/^\#/;
  chomp $line;
  my @F = split (m/\s+/, $line);
  $annot_hash{$F[0]}{hU} = sprintf ("%.1f", $F[3]);
  $annot_hash{$F[0]}{AI} = sprintf ("%.1f", $F[6]); ## round to 1dp
  $annot_hash{$F[0]}{category} = $F[9]; ## key= geneid; val=INGROUP or OUTGROUP
  $annot_hash{$F[0]}{CHS} = $F[10];
  $annot_hash{$F[0]}{tax} = $F[11];
}
close $ANNOT;
print STDERR "[INFO] Collected annotations for ".commify(scalar(keys %annot_hash))." genes\n";

## open groups file
open (my $GROUPS, $orthogroupsfile) or die $!;
while (my $line = <$GROUPS>) {
  chomp $line;
  my @a = split (m/\s+/, $line);
  my @b;
  my $count = 0;
  foreach my $element (@a) {
    ## not all genes in %annot_hash are HGTc...
    if ($annot_hash{$element}) {
      ## test if HGTc; must pass hU > 30 && category eq OUTGROUP with CHS support > 90%
      if ( ($annot_hash{$element}{hU} > 30) and ($annot_hash{$element}{category} eq "OUTGROUP") and ($annot_hash{$element}{CHS} > 90) ) {
        ## is HGTc...
        if ($simple) {
          my $new_id = join (":", $element, $annot_hash{$element}{category}); ## here 'category' should alwys be OUTGROUP
          push (@b, $new_id);
        } else {
          my $new_id = join (":", $element, $annot_hash{$element}{hU}, $annot_hash{$element}{AI}, $annot_hash{$element}{category}, $annot_hash{$element}{CHS}, $annot_hash{$element}{tax});
          push (@b, $new_id);
        }
        $count++;
      } else {
        ## is not HGTc...
        push (@b, $element);
      }
    } else { ## and not all genes are in %annot_hash (eg if there was no hit to UniRef90)
      push (@b, $element);
    }
  }

  ## print
  print $OUT join (" ", @b) . "\n";
  ## (scalar(@b)-1) because the OG id is also an element of @b, so need to discount
  print $PROP join (" ", $b[0], (scalar(@b)-1), $count, (sprintf("%.3f",($count/(scalar(@b)-1))))) . "\n";
}
close $GROUPS;
close $OUT;
close $PROP;

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
