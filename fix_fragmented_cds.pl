#!/usr/bin/env perl

## reubwn December 2022

use strict;
use warnings;
use Getopt::Long;

use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use File::Basename;
use Sort::Naturally;
use List::Util 'sum';
use Data::Dumper;

my $usage = "
SYNOPSIS
  Fix fragmented CDS from poor quality genome annotations using inference from orthology.

OPTIONS [*required]
  -a|--aa          *[FILE]            : target aa sequences (fasta)
  -d|--db          *[DIR]             : path to dir of database proteome files (fasta)
  -g|--orthogroups *[FILE]            : OrthoGroups.txt file (from OrthoFinder)
  -i|--ignore       [STRING[,STRING]] : ID of any taxa to ignore in OF output
  -n|--dna          [FILE]            : target DNA sequences (fasta)
  -o|--out          [STR]             : outfile suffix
  -l|--logfile                        : print stats to logfile [no]
  -h|--help                           : print this message
\n";

my ($aa_file, $dna_file, $db_dir_path, $orthogroups_file, $target_id, $msa_dir, $ignore_string, $logfile, $help);
my $outsuffix = "fixed";

GetOptions (
  'a|aa=s'      => \$aa_file,
  'd|db=s'      => \$db_dir_path,
  'g|orthogroups=s' => \$orthogroups_file,
  't|target_id=s' => \$target_id,
  'm|msa=s'     => \$msa_dir,
  'i|ignore:s'  => \$ignore_string,
  'n|dna:s'     => \$dna_file,
  'o|out:s'     => \$outsuffix,
  'l|logfile'   => \$logfile,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ($target_id && $msa_dir);

my @ignore;
if ( $ignore_string ) {
  @ignore = split (",", $ignore_string);
}

open (my $LOG, ">filter.log") or die $!;
open (my $OUT, ">filter_res.out") or die $!;

my %target_lengths;
my $aa_in = Bio::SeqIO -> new ( -file => $aa_file, -format => 'fasta' );
while (my $seq = $aa_in->next_seq) {
  $target_lengths{$seq->display_id} = $seq->length;
}
print STDERR "[INFO] Number of proteins in '$aa_file': ".scalar(keys %target_lengths)."\n";

my @msa_files = glob ("$msa_dir/*fa");
print STDERR "[INFO] Number of \*.fa MSA files in '$msa_dir' to analyse: ".scalar(@msa_files)."\n";

my $n;
my $m;

foreach my $msa (@msa_files) {
  my $aln_in = Bio::AlignIO -> new ( -file => $msa, -format => 'fasta' );
  my $aln = $aln_in -> next_aln();
  my %counts;
  # my %lengths;
  # my %pairwise_matches;
  my %target_members;
  my $longest_target_length = 0;
  my $longest_target_id;

  foreach my $seq1 ($aln -> each_seq()) {
    my @a = split(/\|/, $seq1->display_id());
    $counts{$a[0]}++;

    # foreach my $seq2 ($aln -> each_seq()) {
    #   if ($seq1->display_id() eq $seq2->display_id()) {
    #     next;
    #   } else {
    #     my $num_matches = ( lc $seq1->seq() ^ lc $seq2->seq() ) =~ tr/\0//;
    #   }
    # }

    if ($a[0] eq $target_id) {
      # print STDERR "[INFO] ".$seq1->display_id.": ".$target_lengths{$seq1->display_id}."\n";
      $target_members{$seq1->display_id} = $target_lengths{$seq1->display_id};
      if ($target_lengths{$seq1->display_id} > $longest_target_length) {
        $longest_target_length = $target_lengths{$seq1->display_id};
        $longest_target_id = $seq1->display_id;
      }
    }
  }

  my %counts_copy = %counts;
  delete $counts_copy{$target_id};
  foreach (@ignore) {
    delete $counts_copy{$_};
  }

  ## number of keys == sum of values, then 1-1 orthogroup (ignoring @ignore)
  if ( ($counts{$target_id}) and (scalar keys %counts_copy >0) ) {
    if ( ($counts{$target_id} > 1) and (scalar keys %counts_copy == sum values %counts_copy) ) {
      print $LOG "[INFO] $msa: $target_id has $counts{$target_id} copies\n";
      print $LOG "[INFO] Longest target in aln is $longest_target_id ($longest_target_length aa)\n";
      delete $target_members{$longest_target_id};
      print $LOG "[INFO] Target PIDs to be removed:\n       ".join(", ", nsort keys %target_members)."\n\n";
      print $OUT join("\n", nsort keys %target_members);
      $m+= scalar(keys %target_members);
      $n++;

    }
  }
}
close $LOG;
close $OUT;

print STDERR "[INFO] Total number of MSAs with >1 copy for $target_id: $n\n";
print STDERR "[INFO] Total number of proteins to be suppressed: $m\n".`date`."\n";

#############

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
    return $rounded;
}