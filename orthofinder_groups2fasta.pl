#!/usr/bin/env perl

## author: reubwn October 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Term::ANSIColor;
use Sort::Naturally;
use File::Path 'rmtree';

my $usage = "
SYNOPSIS
  Converts Orthogroups.txt to fasta.

OPTIONS
  -i|--groups [FILE] : Orthogroups.txt file [required]
  -p|--path   [FILE] : path to directory of fasta files used to construct Orthogroups.txt [required]
  -d|--outdir [STR]  : output dir name [default: orthofinder_groups2fasta]
  -h|--help          : prints this help message
\n";

my ($orthogroups,$path,$outdir,$help);
my $prefix = "orthofinder_groups2fasta";

GetOptions (
  'i|in=s'     => \$orthogroups,
  'p|path=s'   => \$path,
  'd|outdir:s' => \$outdir,
  'h|help'     => \$help,
);

die $usage if $help;
die $usage unless ($orthogroups && $path);

## get sequences
my %seq_hash;
my @fastas = glob("$path*fasta");
if (scalar(@fastas) == 0) {
  print STDERR "[INFO] Nothing found in $fastas with *.fasta... will try *.faa\n";
  @fastas = glob("$fastas*faa");
}
if (scalar(@fastas) == 0) {
  print STDERR "[INFO] Nothing found in $fastas with *.faa... will try *.fna\n";
  @fastas = glob("$fastas*fna");
}
if (scalar(@fastas) == 0) {
  print STDERR "[INFO] Nothing found in $fastas with *.fna... will try *.aa (augustus style)\n";
  @fastas = glob("$fastas*aa");
}
if (scalar(@fastas) == 0) {
  die "[ERROR] Still nothing found in $fastas\nPlease make sure there are protein fastas in $fastas with extension fasta|faa|fna|aa\n";
}
print STDERR "[INFO] Reading sequences from:\n";
foreach (@fastas){
  print STDERR colored($_, 'white on_blue') . "\n";
  my $in = Bio::SeqIO->new ( -file => $_, -format => "fasta" );
  while ( my $seq_obj = $in->next_seq() ){
    $seq_hash{($seq_obj->display_id())} = ($seq_obj->seq());
  }
}
print STDERR "[INFO] Read in ".commify(scalar(keys %seq_hash))." sequences from ".commify(scalar(@fastas))." files\n\n";

## make $outdir
if (-e $outdir && -d $outdir) {
  rmtree([ "$outdir" ]);
  mkdir $outdir;
} else {
  mkdir $outdir;
}

## open Orthogroups
my %orthogroups;
open (my $GROUPS, $orthogroups) or die $!;
while (my $line = <$GROUPS>) {
  chomp $line;
  my @a = split(/\:\s+/, $line);
  my @b = split(/\s+/, $a[1]);
  ## open file:
  open (my $OUT, ">$outdir/$a[0]") or die $!;
  foreach (@b) {
    print $OUT ">$_\n$seq_hash{$_}\n";
  }
  close $OUT;
  print STDERR "\r[INFO] Working on OG: $a[0]"; $|=1;
}

print STDERR "\n[INFO] Finished on ".`date`."\n";
