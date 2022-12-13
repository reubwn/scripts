#!/usr/bin/env perl

## reubwn December 2022

use strict;
use warnings;
use Getopt::Long;

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Fix fragmented CDS from poor quality genome annotations using inference from orthology.

OPTIONS [*required]
  -a|--aa          *[FILE]            : target aa sequences (fasta)
  -d|--db          *[DIR]             : dir of database proteome files (fasta)
  -g|--orthogroups *[FILE]            : OrthoGroups.txt file (from OrthoFinder)
  -i|--ignore       [STRING[,STRING]] : ID of any taxa to ignore in OF output
  -n|--dna          [FILE]            : target DNA sequences (fasta)
  -o|--out          [STR]             : outfile suffix
  -l|--logfile                        : print stats to logfile [no]
  -h|--help                           : print this message
\n";

my ($aa_file, $dna_file, $db_dir, $orthogroups_file, $ignore_string, $logfile, $help);
my $outsuffix = "fixed";

GetOptions (
  'a|aa=s'      => \$aa_file,
  'd|db=s'      => \$db_dir,
  'g|orthogroups=s' => \$orthogroups_file,
  'i|ignore:s'  => \$ignore_string,
  'n|dna:s'     => \$dna_file,
  'o|out:s'     => \$outsuffix,
  'l|logfile'   => \$logfile,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $aa_file && $orthogroups_file && $db_dir );

my @ignore;
if ( $ignore_string ) {
  @ignore = split (",", $ignore_string);
}

my %target_hash;
my $prots_fh = Bio::SeqIO -> new ( -file => $aa_file, -format => 'fasta');
while (my $seq_obj = $prots_fh -> next_seq) {
  $target_hash{$seq_obj->display_id()} = $seq_obj->seq();
}
print STDOUT "[INFO] Number of seqs in '$aa_file': ".commify(keys %target_hash)."\n";


## read in proteins file to hash
## parse Orthogroups.txt file
## find 1-1's, ignoring any GID in 'ignore'
## iterate thru 1-1's, find any with multiple proteins from target
## iterate thru these;

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
