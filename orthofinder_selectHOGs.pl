#!/usr/bin/env perl

## author: reubwn Feb 2022

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS:
  Parses Orthofinder output file 'Orthogroups.GeneCount.tsv' for groups corresponding to input criteria.
  Construct the 'select' string like so: 'ID1=X,ID2=Y,IDN=Z', where ID1 is the unique species ID
  corresponding to the column headers in 'Orthogroups.GeneCount.tsv', and X is the number of members from
  species 1 you want to select groups containing.

OUTPUTS:


OPTIONS:
  -i|--in       [FILE] : Orthogroups.GeneCount.tsv file [required]
  -s|--select   [STR]  : select string, eg. 'ID1=1,ID2=1,IDN=1' would return 1:1 orthologous groups [required]
  -m|--msa_dir  [DIR]  : path to 'MultipleSequenceAlignments' directory of fasta alignments
  -t|--tree_dir [DIR]  : path to 'Resolved_Gene_Trees' directory of tree files
  -o|--out      [STR]  : output dirname [default: 'selected']
  -h|--help            : prints this help message
\n\n";

my ($in,$select,$msa_dir,$tree_dir,$out,$help);
my $outdir = "selected";
my $delim = "|";

GetOptions (
  'i|in=s'     => \$in,
  's|select=s' => \$select,
  'm|msa_dir:s'    => \$msa_dir,
  't|tree_dir:s' => \$tree_dir,
  'o|out:s'    => \$outdir,
  'h|help'     => \$help,
);

die $usage if $help;
die $usage unless ($in && $select);

## parse select string
my %select_hash;
my @a = split(",",$select);
foreach (@a){
  my @b = split("=",$_);
  $select_hash{$b[0]} = $b[1]; ## key= SPECIES_ID; value= PER_SPECIES_GROUP_SIZE
}

mkdir $outdir;

print STDERR "[INFO] Select per-species group sizes:\n";
open (my $INFO, ">$outdir/info.txt") or die $!;
foreach (nsort keys %select_hash){
  print STDERR "- $_ => $select_hash{$_}\n";
  print $INFO "$_\t $select_hash{$_}\n";
}
print STDERR "\n";
close $INFO.

my %table;
my @header;

## open Orthogroups.GeneCount.tsv
open (my $TSV, $in) or die;
while (my $line = <$TSV>) {
  chomp $line;
  my @F = split (/\s+/, $line);
  if ($. == 1) {
    pop @F;
    shift @F;
    @header = @F;
    # print STDERR "@header\n";
    # @table{@F} = (); ## set species IDs as keys
  } else {
    for my $i (0 .. $#header) {
      $table{$F[0]}{$header[$i]} = $F[$i+1] ## key=OG; val={key=species; val=count}
    }
  }
}

my %results;

# print Dumper(\%table);

foreach my $og (nsort keys %table) {
  print STDERR "\r[INFO] $og"; $|=1;
  my $flags = 0;
  foreach my $species (nsort keys %{$table{$og}}) {
    foreach my $select (nsort keys %select_hash) {
      if ( ($select eq $species) && ($table{$og}{$species} == $select_hash{$select}) ) {
        $flags++;
      }
    }
    if ($flags == scalar(keys %select_hash)) {
      $results{$og}++;
    }
  }
}
print STDERR "\n";

print STDERR "[INFO] Found ".scalar(keys %results)." OGs that satisfy selection\n";

open (my $RESULT, ">$outdir/list.txt") or die $!;
foreach (nsort keys %results) {
  print $RESULT "$_\n";

  if ($tree_dir) {
    my $treefile = $_."_tree.txt";
    `cat $tree_dir/$treefile >> $outdir/trees.newick`;
    `printf "\n" >> $outdir/trees.newick`;
  }
}

close $RESULT;
