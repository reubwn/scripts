#!/usr/bin/env perl

## author: reubwn May 2019

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
  -t|--tree       [STRING] : also produce a tree? ('iqtree' or 'fasttree') [default: none]
  -e|--evalue     [STRING] : hmmsearch reporting evalue threshold [1e-5]
     --par_search [STRING] : additional params to be passed to hmmsearch, e.g. --par_search='--domE 1e-20'
     --par_select [STRING] : additional regexp to be passed to fastaqual_select.pl, e.g. --par_select='-regexp t1'
  -x|--excl_alt            : exclude alternative transcripts from the same locus
  -h|--help                : prints this help message
\n";

my ($query,$db,$ali,$tree,$exclude_alt_same_locus,$help);
my $par_search = "";
my $par_select = "";
my $evalue = "1e-5";
my $prefix = "hmmsearch-easy";

GetOptions (
  'q|query=s'     => \$query,
  'd|db=s'        => \$db,
  'a|ali=s'       => \$ali,
  't|tree:s'      => \$tree,
  'e|evalue:s'    => \$evalue,
  'par_search:s'  => \$par_search,
  'par_select:s'  => \$par_select,
  'x|excl_alt'    => \$exclude_alt_same_locus,
  'p|prefix:s'    => \$prefix,
  'h|help'        => \$help,
);

die $usage if $help;
die $usage unless ($query && $db && $ali);

if (system("hmmsearch -h &> /dev/null")!=0) {
  die "[ERROR] hmmsearch returned non-zero exit status!\n";
  die "[ERROR] Is hmmsearch installed and in \$PATH?\n";
}

if (system("fastaqual_select.pl &>/dev/null")!=0) {
  die "[ERROR] fastaqual_select.pl returned non-zero exit status!\n";
  die "[ERROR] Is fastaqual_select.pl installed and in \$PATH?\n";
}

die "[ERROR] File $query, $db or $ali does not exist!\n" if ((! -f $query) or (! -f $db) or (! -f $ali));

print STDERR "[INFO] Query alignment: $ali\n";
print STDERR "[INFO] Query HMM file: $query\n";
print STDERR "[INFO] Subject database: $db\n";
print STDERR "[INFO] Additional hmmsearch parameters: '$par_search'\n" if $par_search =~ m/\w+/;
print STDERR "[INFO] Additional fastaqual_select.pl parameters: '$par_select'\n" if $par_select =~ m/\w+/;

## run hmmsearch
`hmmsearch -o $prefix.hmmsearch.$evalue.out --tblout $prefix.hmmsearch.$evalue.tblout --noali -E $evalue $par_search $query $db`;

## check results and die if none
my $result = `grep -v "^#" $prefix.hmmsearch.$evalue.tblout`;
unless ($result) {
  die "[INFO] No matches found! Stopping here\n";
}

## get hits to faa
`awk '{print \$1}' $prefix.hmmsearch.$evalue.tblout | fastaqual_select.pl -f $db -i - $par_select > $prefix.hmmsearch.$evalue.out.faa`;

## exclude alternative transcripts from the same locus
if ($exclude_alt_same_locus) {
  my $in = Bio::SeqIO -> new ( -file => "$prefix.hmmsearch.$evalue.out.faa", -format => "fasta" );
  my $out = Bio::SeqIO -> new ( -file => ">TMP.faa", -format => "fasta" );
  my (%seen, %keep);
  while ( my $seq_obj = $in->next_seq() ) {
    ## assumes BRAKER style alt transcript names convention!
    my $name = $seq_obj->display_id(); $name =~ s/\.t\d+//;
    # print STDERR "$name\n";
    if ($seen{$name}) {
      print STDERR "[INFO] Discarded ".$seq_obj->display_id()." because an alternative transcript for this locus already exists\n";
    } else {
      $out->write_seq($seq_obj);
    }
    $seen{$name}++; ## assumes BRAKER style alt transcript names convention!
  }

  `mv TMP.faa $prefix.hmmsearch.$evalue.out.faa`;
}

## run hmmalign
`hmmalign -o $prefix.hmmalign.sto --mapali $ali --trim $query $prefix.hmmsearch.$evalue.out.faa`;

## convert aligned *.sto to clustal
`esl-reformat clustal $prefix.hmmalign.sto > $prefix.hmmalign.clustal`;

## then convert clustal to fasta
`seqret -sequence $prefix.hmmalign.clustal -outseq $prefix.hmmalign.fasta &> /dev/null`;

## make a tree using FastTree or IQ-TREE if desired
if ( $tree ) {
  if ( $tree eq "fasttree" ) {
    if ( system("fasttree &> /dev/null")!=0 ) {
      die "[ERROR] FastTree returned non-zero exit status!\n";
      die "[ERROR] Is FastTree installed and in \$PATH?\n";
    } else {
      `fasttree $prefix.hmmalign.fasta > $prefix.hmmalign.fasta.tree`;
    }

  } elsif ( $tree eq "iqtree" ) {
    if ( system("iqtree -h &> /dev/null")!=0 ) {
      die "[ERROR] IQ-TREE returned non-zero exit status!\n";
      die "[ERROR] Is IQ-TREE installed and in \$PATH?\n";
    } else {
      `iqtree -s $prefix.hmmalign.fasta -nt 4 -bb 1000 -m TEST`;
    }
  } else {
    die "[ERROR] Did not recognise tree command, please choose either 'fasttree' or 'iqtree'\n";
  }
}

print STDERR "[INFO] Done " . `date`;
