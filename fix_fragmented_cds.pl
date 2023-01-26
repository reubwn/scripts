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
# die $usage unless ( $aa_file && $orthogroups_file && $db_dir_path );

# my @prots_files = glob("$db_dir_path/*fa $db_dir_path/*faa $db_dir_path/*fasta");
# my $total_proteomes = scalar(@prots_files);
# print STDERR "[INFO] Path to proteome files: $db_dir_path\n";
# print STDERR "[INFO] Number of input files found there: $total_proteomes\n";

my @ignore;
if ( $ignore_string ) {
  @ignore = split (",", $ignore_string);
}

my @msa_files = glob ("$msa_dir/*fa");
print STDERR "[INFO] Number of \*.fa MSA files in '$msa_dir': ".scalar(@msa_files)."\n";

my $n;
foreach my $msa (@msa_files) {
  my $in = Bio::AlignIO -> new ( -file => $msa, -format => 'fasta' );
  my $aln = $in -> next_aln();
  my %counts;
  my %lengths;
  my %pairwise_matches;
  foreach my $seq1 ($aln -> each_seq()) {
    my @a = split(/\|/, $seq1->display_id());
    $counts{$a[0]}++;
    push (@{$lengths{$a[0]}{pids}}, $seq1->display_id());
    push (@{$lengths{$a[0]}{lengths}}, $seq1->length());
    foreach my $seq2 ($aln -> each_seq()) {
      if ($seq1->display_id() eq $seq2->display_id()) {
        next;
      } else {
        my $num_matches = ( lc $seq1->seq() ^ lc $seq2->seq() ) =~ tr/\0//;
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
      print STDERR "[INFO] $msa: $target_id has $counts{$target_id} copies\n";
      $n++;

    }
  }
}
print STDERR "[INFO] Total number of MSAs with >1 copy for $target_id: $n\n";




#
#
# my %target_hash;
# my $prots_fh = Bio::SeqIO -> new ( -file => $aa_file, -format => 'fasta');
# while (my $seq_obj = $prots_fh -> next_seq) {
#   $target_hash{$seq_obj->display_id()} = $seq_obj->seq();
# }
# print STDERR "[INFO] Number of seqs in '$aa_file': ".commify(scalar(keys %target_hash))."\n";
#
# ## read in proteins file to hash
# ## parse Orthogroups.txt file
# ## find 1-1's, ignoring any GID in 'ignore'
# ## iterate thru 1-1's, find any with multiple proteins from target
# ## iterate thru these;
#
# open (my $ORTHO, $orthogroups_file) or die $!;
# while (my $line = <$ORTHO>) {
#   chomp $line;
#   my @a = split (/:\s/, $line);
#   my @b = split (/\s+/, $a[1]);
#   my %gids;
#   foreach my $entry (@b) {
#     my @c = split (/\|/, $entry);
#     $gids{$c[0]}++;
#   }
#   ## delete the target gids which we don't want to consider here
#   foreach (@ignore) {
#     delete ($gids{$_});
#   }
#   if ( (scalar(keys %gids) == $total_proteomes) and (sum values %gids == $total_proteomes-(scalar(@ignore)))) {
#     print "$line\n";
#   }
# }

## iterate thru all but @ignore and pick out 1-1's

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
