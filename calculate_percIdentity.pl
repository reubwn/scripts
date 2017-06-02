#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Bio::AlignIO;
use Sort::Naturally;
use Bio::SimpleAlign;

my $usage = "
SYNOPSIS
  Calculates average percentage identity for a dir of alignments.

  Treats gaps in 4 ways:
    1. PID: ignores any gap vs base or gap vs gap site
    2. PID_ALN: counts number of identical columns, divides by total aln length (incl gaps)
    3. PID_SHORT: counts number of identical columns, divides by shortest seq length (excl gaps)
    4. PID_LONG: counts number of identical columns, divides by longest seq length (excl gaps, but penalises unaligned columns of longest seq)

  Assumes fasta alignments.

OPTIONS
  -i|--indir [DIR]    : path to dir of alignments
  -o|--out   [STRING] : name of outfile (default = 'perc_ident.txt')
  -h|--help           : this message
\n";

my ($indir, $help);
my $outfile = "perc_ident.txt";

GetOptions (
  'indir|i=s' => \$indir,
  'out|o:s'   => \$outfile,
  'help|h'    => \$help
);

die $usage if $help;
die $usage unless ($indir);

## get files:
my @files=<$indir/*>;
print STDERR "[INFO] Number of files: ".scalar(@files)."\n";

## open out:
open (my $OUT, ">$outfile") or die "[ERROR] Cannot open $outfile: $!\n";
print $OUT "FILE\tPID\tPID_ALN\tPID_SHORT\tPID_LONG\n";

foreach (nsort @files) {
  my $filename = $_;
  $filename =~ s/.+\///;
  print STDERR "[INFO] Processing file: $filename\r"; $|=1;
  my $in = Bio::AlignIO->new ( -file=>$_, -format=>"fasta");
  my $alnobj = $in->next_aln();
  print $OUT "$filename\t".($alnobj->percentage_identity())."\t".($alnobj->overall_percentage_identity('align'))."\t".($alnobj->overall_percentage_identity('short'))."\t".($alnobj->overall_percentage_identity('long'))."\n";
}
print STDERR "\n[INFO] Finished on ".`date`."\n";
