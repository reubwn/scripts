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
  -a|--annot       [FILE] : *HGT_results file
  -o|--outfile     [STR]  : output base filename (default: 'inputfilename')
  -h|--help               : print this message
\n";

my ($orthogroupsfile, $annot, $outfile, $help);

GetOptions (
  'i|orthogroups=s' => \$orthogroupsfile,
  'a|annot=s'       => \$annot,
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
while (my $line = <$ANNOT>) {
  chomp $line;
  my @F = split (m/\s+/, $line);
  $annot_hash{$F[0]}{hU} = $F[3];
  $annot_hash{$F[0]}{AI} = (sprintf "%.1f", $F[6]); ## round to 1dp
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
  foreach my $element (@a) {
    if ($annot_hash{$element}) {
      if ($annot_hash{$element}{category} eq "OUTGROUP") {
        my $new_id = join (":", $element, $annot_hash{$element}{hU}, $annot_hash{$element}{AI}, $annot_hash{$element}{category}, $annot_hash{$element}{CHS}, $annot_hash{$element}{tax});
        push (@b, $new_id);
      } else {
        my $new_id = join (":", $element, $annot_hash{$element}{hU}, $annot_hash{$element}{category});
        push (@b, $new_id);
      }
    } else {
      push (@b, $element);
    }
  }
  ## count number OUTGROUP seqs in OG
  my $count = () = join(" ", @b) =~ m/OUTGROUP/g;

  ## print
  print $OUT join (" ", @b) . "\n";
  print $PROP join (" ", scalar(@b), $count, (sprintf("%.3f",($count/scalar(@b))))) . "\n";
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
