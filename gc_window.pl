#!/usr/bin/env perl

## author: reubwn October 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
SYNOPSIS:

OUTPUTS:

OPTIONS:
  -f|--fasta     [FILE]  : input fasta file [required]
  -w|--window    [INT]   : window size to calculate \%GC over [default: 500]
  -s|--step      [INT]   : step size for next window [default: same as window size]
  -N|--threshold [FLOAT] : threshold for proportion of NNNs in window, otherwise will print \"NA\" [default: 0.5]
  -o|--out       [STR]   : output filename [default: print to STDOUT]
  -n|--noheader          : omit header from output [default: print header]
  -h|--help              : prints this help message

EXAMPLES:

\n";

my ($fasta,$out,$noheader,$help);
my $window = 500;
my $step = $window;
my $threshold_NNNs = 0.5;

GetOptions (
  'fastas|f=s'    => \$fasta,
  'window|w:s'    => \$window,
  'step|s:s'      => \$step,
  'threshold|N:f' => \$threshold_NNNs,
  'out|o:s'       => \$out,
  'noheader|n'    => \$noheader,
  'help|h'        => \$help,
);

die $usage if $help;
die $usage unless $fasta;

## print to file if specified otherwise to STDOUT
my $OUT;
if ($out) {
  open ($OUT, ">$out") or die "$!\n";
} else {
  $OUT = \*STDOUT;
}
print $OUT join "\t", "seq_name", "window", "start", "end", "window", "num_gc", "num_nnn", "prop_gc", "\n" unless $noheader;

my %GC;
my ($start,$end) = (1,$window); ## seqstrings are 1-offset, ie first base is at pos 1

my $in = Bio::SeqIO -> new( -file => $fasta, -format => "fasta" );
while (my $seq = $in -> next_seq()){
  my $len = $seq -> length();
  my $w = 1;
  for (my $i=$window; $i<=$len; $i+=$step){
    my $substring = $seq -> subseq($start,$window);
    my $gc = $substring =~ tr/GCgc//; ## count GC
    my $nn = $substring =~ tr/Nn//; ## count NNNs
    my ($prop_gc, $prop_nn);
    if ((length($substring)-$nn) == 0) { ## window is all NNNs
      $prop_gc = "NA";
      $prop_nn = 1;
    } elsif ($nn > 0) { ## there are some NNNs in window
      $prop_nn = sprintf("%.4f",($nn/$window);
      if ($prop_nn >= $threshold_NNNs) { ## if proportion NNNs is greater than threshold
        $prop_gc = "NA";
      } else {
        $prop_gc = sprintf("%.4f",($gc/(length($substring)-$nn)));
      }
    } else { ## there are no NNNs
      $prop_gc = sprintf("%.4f",($gc/$window);
    }
    print $OUT join "\t", $seq->display_id, $w, commify($start), commify($window), commify(length($substring)), commify($gc), commify($nn), $prop_gc, "\n";
    $start += $step;
    $window += $step;
    $w++;
  }
}

############################# subs

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
