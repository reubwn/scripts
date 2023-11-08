#!/usr/bin/env perl

## reubwn Nov 23

use strict;
use warnings;
use Getopt::Long;

use File::Path qw( rmtree );
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Collapse per-gene CPM values based on orthogroups.

OPTIONS [*required]
  -i|--in          *[FILE] : FOFN of files to parse column names/values from
  -g|--orthogroups *[FILE] : Orthogroups.txt file
  -s|--species       [STR] : species IDs to pull out of Orthogroups file (comma delim)
  -o|--out           [STR] : outfiles prefix ('result')
  -h|--help                : print this message
\n";

my ($data_fofn, $orthogroups_file, $help);
my $species_input = "all";
my $outprefix = "result";

GetOptions (
  'i|in=s'      => \$data_fofn,
  'g|orthogroups=s' => \$orthogroups_file,
  's|species:s' => \$species_input,
  'o|out:s'     => \$outprefix,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $data_fofn && $orthogroups_file );

## parse species IDs from input
my @species_ids = split (/\,/, $species_list);
print STDERR "[INFO] Species to select from orthogroups: ".join(", ", @species_ids) . "\n";

## read in OGs
my %orthogroups;
open (my $ORTHO, $orthogroups_file) or die $!;
while (my $line = <$ORTHO>) {
  chomp $line;
  my @F = split (m/\s+/, $line);
  my $og_id = shift @F;
  ## pull out gene ids for individual species based on OG convention e.g. "AVAG|g67330.t1"
  my @arr1 = grep (/^ARIC\|/, @F);
  my @arr2 = grep (/^AVAG\|/, @F);
  ## push into hash
  foreach my $species_id (@species_ids) {
    print STDERR "[INFO] Pulling out genes for '$species_id'...\n";
    $orthogroups{$species_id}{$og_id} = \@F;
  }


}
close $ORTHO;
