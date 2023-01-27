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
  -m|--msa         *[DIR]             : path to dir of MSA files (fasta)
  -t|--target      *[STRING]          : ID of target genome
  -i|--ignore       [STRING[,STRING]] : ID of any taxa to ignore
  -o|--out          [STR]             : outfile prefix ('filter')
  -n|--min          [INT]             : mimimum number of species in OG to be included [5]
  -f|--fuzzy        [INT]             : allow up to INT extra copies in non-target species (i.e. departure from strict 1-1's)
  -h|--help                           : print this message
\n";

my ($aa_file, $msa_dir, $target_id, $ignore_string, $fuzzy, $help);
my $outprefix = "filter";
my $min_OG_size = 5;

GetOptions (
  'a|aa=s'      => \$aa_file,
  't|target=s'  => \$target_id,
  'm|msa=s'     => \$msa_dir,
  'i|ignore:s'  => \$ignore_string,
  'o|out:s'     => \$outprefix,
  'f|fuzzy:i'   => \$fuzzy,
  'n|min:i'     => \$min_OG_size,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $aa_file && $msa_dir && $target_id );

## parse ignore string
my @ignore;
if ( $ignore_string ) {
  @ignore = split (",", $ignore_string);
}

## open out files
open (my $LOG, ">$outprefix.log") or die $!;
open (my $OUT, ">$outprefix"."_exclude.out") or die $!;

## parse aa target file
my %target_lengths;
my $aa_in = Bio::SeqIO -> new ( -file => $aa_file, -format => 'fasta' );
while (my $seq = $aa_in->next_seq) {
  $target_lengths{$seq->display_id} = $seq->length;
}
print STDERR "[INFO] Number of proteins in target file '$aa_file': ".scalar(keys %target_lengths)."\n";

## glob MSA files
my @msa_files = glob ("$msa_dir/*fa");
print STDERR "[INFO] Number of \*.fa MSA files in '$msa_dir' to analyse: ".scalar(@msa_files)."\n";

my $n;
my $m;

## iterate thru MSA files
foreach my $msa (@msa_files) {
  my $aln_in = Bio::AlignIO -> new ( -file => $msa, -format => 'fasta' );
  my $aln = $aln_in -> next_aln();
  my %counts;
  # my %lengths;
  # my %pairwise_matches;
  my %target_members;
  my $longest_target_length = 0;
  my $longest_target_id;

  ## iterate thru seqs in current MSA
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

    ## find the longest target seq in current MSA/OG
    if ($a[0] eq $target_id) {
      # print STDERR "[INFO] ".$seq1->display_id.": ".$target_lengths{$seq1->display_id}."\n";
      $target_members{$seq1->display_id} = $target_lengths{$seq1->display_id};
      if ($target_lengths{$seq1->display_id} > $longest_target_length) {
        $longest_target_length = $target_lengths{$seq1->display_id};
        $longest_target_id = $seq1->display_id;
      }
    }
  }

  ## make a copy of counts
  my %counts_copy = %counts;
  ## remove target counts
  delete $counts_copy{$target_id};
  ## remove any other GIDs you don't want to include
  foreach (@ignore) {
    delete $counts_copy{$_};
  }

  ## skip MSAs that do not contain any seqs from target, or that are composed entirely of seqs from target or @ignore
  if ( ($counts{$target_id}) and (scalar keys %counts_copy > 0) ) {

    if ( $fuzzy ) {
      if ( ($counts{$target_id} > 1) and (scalar keys %counts_copy >= $min_OG_size) and (sum values %counts_copy <= (scalar keys %counts_copy + $fuzzy)) ) {
        print $LOG "[INFO] $msa: $target_id has $counts{$target_id} copies\n";
        print $LOG "[INFO] Longest target in aln is $longest_target_id ($longest_target_length aa)\n";
        delete $target_members{$longest_target_id};
        print $LOG "[INFO] Target PIDs to be removed:\n       ".join(", ", nsort keys %target_members)."\n\n";
        print $OUT join("\n", nsort keys %target_members)."\n";
        $m+= scalar(keys %target_members);
        $n++;
        
      }
    } else {
      ## number of keys == sum of values, then 1-1 orthogroup (ignoring @ignore)
      if ( ($counts{$target_id} > 1) and (scalar keys %counts_copy >= $min_OG_size) and (sum values %counts_copy == scalar keys %counts_copy) ) {
        print $LOG "[INFO] $msa: $target_id has $counts{$target_id} copies\n";
        print $LOG "[INFO] Longest target in aln is $longest_target_id ($longest_target_length aa)\n";
        delete $target_members{$longest_target_id};
        print $LOG "[INFO] Target PIDs to be removed:\n       ".join(", ", nsort keys %target_members)."\n\n";
        print $OUT join("\n", nsort keys %target_members)."\n";
        $m+= scalar(keys %target_members);
        $n++;

      }
    }
  }
}
close $LOG;
close $OUT;

print STDERR "[INFO] Total number of MSAs with >1 copy for $target_id: $n\n";
print STDERR "[INFO] Total number of proteins to be suppressed: $m\n[INFO] Finished ".`date`."\n";

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
