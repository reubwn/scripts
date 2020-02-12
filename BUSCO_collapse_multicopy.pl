#!/usr/bin/env perl

## author: reubwn Feb 2020

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use Data::Dumper;

my $usage = "
SYNOPSIS:
  Collapse multiple BUSCO copies to single copy by random selection

OPTIONS:
  -b|--busco [PATH] : path to BUSCO results directory [required]
  -h|--help         : prints this help message
\n";

my ($busco_path, $help, $debug);

GetOptions (
  'b|busco=s' => \$busco_path,
  'h|help'    => \$help,
  'debug'     => \$debug
);

die $usage if $help;
die $usage unless ($busco_path);

############
## get prots
############

my %extracted_proteins_hash;
my @faa_files = glob ("$busco_path/augustus_output/extracted_proteins/*faa*");
foreach my $faa_file (@faa_files) {
  print STDERR "\r[INFO] Extracting proteins from $faa_file"; $|=1;
  my $in = Bio::SeqIO -> new ( -file => $faa_file, -format => "fasta");
  while ( my $seq_obj = $in -> next_seq() ) {
    $extracted_proteins_hash{$seq_obj->display_id()} = $seq_obj->seq(); ##key= seq name; val= seq
  }
}
print STDERR "\n[INFO] Extracted ".scalar(keys %extracted_proteins_hash)." proteins from $busco_path/augustus_output/extracted_proteins/\n";

####################
## parse BUSCO table
####################

my @full_table_file = glob ("$busco_path/full_table*tsv");
if (@full_table_file != 1) {
  die "[ERROR] There are either zero or multiple full_table files!\n";
}

my $count = 1;
my %full_table_hash;
open (my $FULL, $full_table_file[0]) or die $!;
while (<$FULL>) {
  if (m/^\#/) {
    next;
  } else {
    chomp;
    my @F = split (m/\s+/, $_);
    my %busco_hash = (
      $F[0] => {
        'Status' => $F[1],
        'Contig' => $F[2],
        'Start' => $F[3],
        'End' => $F[4],
        'Score' => $F[5],
        'Length' => $F[6]
      }
    );
    push ( @{ $full_table_hash{$F[0]}{Name} }, $F[0] );
    push ( @{ $full_table_hash{$F[0]}{Status} }, $F[1] );
    push ( @{ $full_table_hash{$F[0]}{Contig} }, $F[2] );
    push ( @{ $full_table_hash{$F[0]}{Start} }, $F[3] );
    push ( @{ $full_table_hash{$F[0]}{End} }, $F[4] );
    push ( @{ $full_table_hash{$F[0]}{Score} }, $F[5] );
    push ( @{ $full_table_hash{$F[0]}{Length} }, $F[6] );

  }
}

print STDERR "[INFO] Number keys in hash: ".scalar(keys %full_table_hash)."\n";
## Dumper
print Dumper(\%full_table_hash) if $debug;
