#!/usr/bin/env perl

## author: reubwn May 2020

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;

my $usage = "
SYNOPSIS:
  Runs hmmsearch vs fasta database, and returns alignment (clustal) including
  seeds in the HMM file.

OPTIONS:
  -p|--prefix     [STRING] : output file prefix ['hmmsearch-easy']
  -q|--query      [FILE]   : query file (HMM format) [required]
  -d|--db         [FILE]   : sequence database (fasta format) [required]
  -a|--ali        [FILE]   : query file (aligned .sto or .fasta) [required]
  -e|--evalue     [STRING] : hmmsearch reporting evalue threshold [1e-5]
  -c|--cpu        [INT]    : number of worker threads for ClustalO/IQ-TREE [4]
  -t|--tree                : also produce an IQ-TREE phylogeny?
     --par_search [STRING] : additional params to be passed to hmmsearch
     --par_select [STRING] : additional regexp to be passed to fastaqual_select.pl
  -x|--excl_alt            : only include BRAKER-style primary transcripts ['t1']
  -h|--help                : prints this help message
\n";

my ($query_string,$db,$exclude_alt_same_locus,$ali,$tree,$help);
my $par_search = "";
my $par_select = "";
my $evalue = "1e-5";
my $threads = 8;
my $prefix = "hmmsearch-easy";

GetOptions (
  'q|query=s'     => \$query_string,
  'd|db=s'        => \$db,
  'a|ali:s'       => \$ali,
  't|tree'        => \$tree,
  'e|evalue:s'    => \$evalue,
  'c|cpu:i'       => \$threads,
  'par_search:s'  => \$par_search,
  'par_select:s'  => \$par_select,
  'x|excl_alt'    => \$exclude_alt_same_locus,
  'p|prefix:s'    => \$prefix,
  'h|help'        => \$help,
);

die $usage if $help;
die $usage unless ($query_string && $db);

if (system("hmmsearch -h &> /dev/null")!=0) {
  die "[ERROR] hmmsearch returned non-zero exit status!\n";
  die "[ERROR] Is hmmsearch installed and in \$PATH?\n";
}

if (system("fastaqual_select.pl &>/dev/null")!=0) {
  die "[ERROR] fastaqual_select.pl returned non-zero exit status!\n";
  die "[ERROR] Is fastaqual_select.pl installed and in \$PATH?\n";
}

my @queries = split (m/\,/, $query_string); ## multi-query

die "[ERROR] HMM file $db does not exist!\n" if (! -f $db);

print STDERR "[INFO] Query HMM file(s): @queries\n";
print STDERR "[INFO] Subject database: $db\n";
print STDERR "[INFO] Additional hmmsearch parameters: '$par_search'\n" if $par_search =~ m/\w+/;
print STDERR "[INFO] Additional fastaqual_select.pl parameters: '$par_select'\n" if $par_select =~ m/\w+/;
my @info;

## run hmmsearch over each HMM in turn
foreach my $query (@queries) {
  print STDERR "[INFO] HMM $query vs datatbase $db\n";
  push (@info, $query);
  my $string = join (".", @info);
  $string =~ s/.hmm//g;

  ## exclude BRAKER-style alt transcripts
  if ($exclude_alt_same_locus) {
    print STDERR "[INFO] Including only 't1' proteins (temp.faa)\n";
    `fastaqual_select.pl -f $db -regexp t1 $par_select > temp_db.faa`;
    `hmmsearch -o $prefix.hmmsearch.$evalue.$string.out --tblout $prefix.hmmsearch.$evalue.$string.tblout --noali -E $evalue $par_search $query temp_db.faa`;
  } else {
    `hmmsearch -o $prefix.hmmsearch.$evalue.$string.out --tblout $prefix.hmmsearch.$evalue.$string.tblout --noali -E $evalue $par_search $query $db`;
  }

  ## check results and die if none
  my $wc = `grep -v "^#" $prefix.hmmsearch.$evalue.$string.tblout | wc -l`;
  my $result = `grep -v "^#" $prefix.hmmsearch.$evalue.$string.tblout`;
  chomp $wc;
  chomp $result;
  unless ($result) {
    die "[RESULT] No matches found! Stopping here\n";
  } else {
    print STDERR "[RESULT] Found $wc matches with domains ".join(" and ",@info) . "\n";
    print STDERR "$result\n";
  }

  ## get hits to faa
  `awk '{print \$1}' $prefix.hmmsearch.$evalue.$string.tblout | fastaqual_select.pl -f $db -i - $par_select > $prefix.hmmsearch.$evalue.$string.out.faa`;
  
  ## result inherits $db
  $db = "$prefix.hmmsearch.$evalue.$string.out.faa";
  
}
## delete temp
if ( $exclude_alt_same_locus ) {
  `rm temp_db.faa`;
}

## align final out with clustalo
if ( $ali ) {
  `cat $ali $db > temp_ali.faa`;
  `clustalo --threads=$threads --force -i temp_ali.faa -o $prefix.hmmsearch.$evalue.clustalo.fasta`;
  `rm temp_ali.faa`;
}

## build tree
if ( $ali && $tree ) {
  print STDERR "[INFO] Running IQ-TREE on $prefix.hmmsearch.$evalue.clustalo.fasta...\n";
  `iqtree -s $prefix.hmmsearch.$evalue.clustalo.fasta -nt AUTO -ntmax $threads -bb 1000 -m TEST -msub nuclear -quiet`;
}

print STDERR "[INFO] Done " . `date`;
