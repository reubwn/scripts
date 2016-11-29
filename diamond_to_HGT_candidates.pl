#!/usr/bin/env perl

## author: reubwn Nov 2016

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use List::Util qw(reduce);

my $usage = "
SYNOPSIS:
  Goal is to take a taxified Diamond BLAST file, and for each hit recurse up the
  tax tree until that hit can be categorised into \"ingroup\" versus \"outgroup\"
  (e.g., \"metazoan\" or \"non-metazoan\" etc.) Sum the bitscores across all hits to calculate
  the winner, but only accept the winner if there is a high degree of congruence
  across all hits (e.g., >=90% of all hits agree on the winner).

OUTPUTS:

OPTIONS:
  -i|--in              [FILE] : tab formatted Diamond output file [required]
  -n|--nodesDB         [FILE] : nodesDB.txt file from blobtools [required]
  -t|--taxid_threshold [INT]  : NCBI taxid to recurse up to; i.e., threshold taxid to define \"ingroup\" [default = 33208 (metazoa)]
  -c|--taxid_column    [INT]  : define taxid column for --in (first column = 1) [default: 13]
  -b|--bitscore_column [INT]  : define bitscore column for --in (first column = 1) [default: 12]
  -p|--prefix          [FILE] : filename prefix for outfile [default = INFILE.HGT_decisions.txt]
  -#|--header                 : don't print header [default: do print it]
  -v|--verbose                : say more things [default: be quiet]
  -h|--help                   : prints this help message

EXAMPLES:

\n";

my ($in,$nodesDB,$prefix,$outfile,$header,$verbose,$help);
my $tax_threshold = 33208;
my $tax_column = 13;
my $bitscore_column = 12;
#my $proportion = 0.9;

GetOptions (
  'in|i=s'              => \$in,
  'nodesDB|n=s'         => \$nodesDB,
  'tax_threshold|t:i'   => \$tax_threshold,
  'tax_column|c:i'      => \$tax_column,
  'bitscore_column|b:i' => \$bitscore_column,
  'prefix|p:s'          => \$prefix,
  'header|#'            => \$header,
  'verbose|v'           => \$verbose,
  'help|h'              => \$help,
);

die $usage if $help;
die $usage unless ($in && $nodesDB);

## define outfile:
if ($prefix) {
  $outfile = $prefix.".HGT_decisions.txt";
} else {
  $outfile = $in.".HGT_decisions.txt";
}

## parse nodesDB:
print STDOUT "Building taxonomy databases from $nodesDB...\n";
my (%nodes_hash, %names_hash, %rank_hash, %file_hash);
open(my $NODES, $nodesDB) or die $!;
while (<$NODES>) {
  chomp;
  next if /\#/;
  my @F = split (/\t/, $_);
  $nodes_hash{$F[0]} = $F[3]; ## key= child taxid; value= parent taxid
  $names_hash{$F[0]} = $F[2]; ## key= taxid; value= species name
  $rank_hash{$F[0]} = $F[1]; ## key= taxid; value= rank
}
close $NODES;
print STDOUT "  Nodes parsed: ".scalar(keys %nodes_hash)."\n";
#foreach (nsort keys %tax_hash) { print "$_\t$tax_hash{$_}\n"; }

## test recursive walk:
#print "Testing...\n";
#my $test=317; ## some taxid
#print "Result is: ".tax_walk_to_get_rank($test)."\n";

## parse Diamond file:
print STDOUT "Parsing Diamond file: $in...\n";
my (%bitscore_hash,%sum_bitscores_per_query_hash,%num_hits_per_query_hash);
open (my $DIAMOND, $in) or die $!;
while (<$DIAMOND>) {
  chomp;
  next if /^\#/;
  my @F = split (/\s+/, $_);
  if ($F[($tax_column-1)] !~ m/\d+/) {
    print STDERR "[WARN] The taxid ".$F[($tax_column-1)]." on line $. of $in does not look like a valid NCBI taxid... Skipping this entry\n";
    next;
  }
  $bitscore_hash{$F[($tax_column-1)]} += $F[($bitscore_column-1)]; ## sum bitscore per taxid; key= taxid, value= sumofbitscores
  $sum_bitscores_per_query_hash{$F[0]} = \%bitscore_hash; ## key= query name; value= hash of {key= taxid; value= sum of bitscores per taxid}
  $num_hits_per_query_hash{$F[0]}++; ## not sure if this is needed?
}
close $DIAMOND;
print STDOUT "  Done\n";

## open outfile:
open (my $OUT, ">",$outfile) or die $!;
print $OUT "\#query\ttaxid\tbestsumbitscore\tsuperkingdom;kingdom;phylum\tingrouptaxname\tdecision\tsupport\n" unless $header;

############################################ MAIN CODE

## get winning bitscore and taxid; calculate congruence among all taxids for all hits per query:
foreach my $query (nsort keys %sum_bitscores_per_query_hash) {
  my %bitscore_hash = %{ $sum_bitscores_per_query_hash{$query} }; ## key= taxid; value= summed bitscore
  my (%count_categories, %support_categories);

  ## get support for winning taxid from other hits:
  foreach my $taxid (keys %bitscore_hash) {
    $count_categories{tax_walk($taxid)}++; ## count categories; if each hit's taxid falls within/outwith the $tax_threshold
    print join "\t", $query, $bitscore_hash{$taxid}, $taxid, $names_hash{$taxid}, tax_walk($taxid), "\n" if $verbose;
  }
  foreach my $cat (keys %count_categories) {
    $support_categories{$cat} = percentage($count_categories{$cat}, scalar(keys %bitscore_hash)); ## calculate proportion of support for the category of the winner
    print join "\t", $cat, percentage($count_categories{$cat}, scalar(keys %bitscore_hash))."\%", "\n" if $verbose;
  }

  ## get taxid with highest bitscore:
  my $taxid_with_highest_bitscore = List::Util::reduce { $bitscore_hash{$b} > $bitscore_hash{$a} ? $b : $a } keys %bitscore_hash; ## winning taxid
  my $taxid_with_highest_bitscore_category = tax_walk($taxid_with_highest_bitscore); ## category of winning taxid ("ingroup", "outgroup" or "unassigned")
  my $taxid_with_highest_bitscore_category_support = $support_categories{$taxid_with_highest_bitscore_category}; ## % support from other hits

  ## print to $out
  print $OUT join "\t", $query, $taxid_with_highest_bitscore, $bitscore_hash{$taxid_with_highest_bitscore}, tax_walk_to_get_rank($taxid_with_highest_bitscore), "ingroup=".$names_hash{$tax_threshold}, $taxid_with_highest_bitscore_category, $taxid_with_highest_bitscore_category_support, "\n";
}
close $OUT;

############################################ SUBS

sub tax_walk {
    my $taxid = $_[0];
    my $walk_to;
    if (exists $_[1]) {
      $walk_to = $_[1];
    } else {
      $walk_to = $tax_threshold; ## default is metazoa
    }

    ## first parent:
    my $parent = $nodes_hash{$taxid};

    ## return "unassigned" if hit has no valid taxid
    my $result;
    if ($parent !~ m/\d+/) {
      $result = "unassigned";
      return $result;
    }
#    print "First parent is: $parent\n";

    ## recurse the tree:
    while (1) {
      if ($parent == $walk_to) {
        ## is ingroup
        $result = "ingroup";
        last;

      } elsif ($parent == 1) {
        ## root; i.e., the whole tree has been recursed without finding $threshold, therefore $taxid must reside in another part of the tax tree
        $result = "outgroup";
        last;

      } elsif ($parent == 32644) {
        ## taxid for "unidentified"
        $result = "unassigned";
        last;

      } else {
        ## walk up the tree!
        $parent = $nodes_hash{$parent};
#        print "  Parent is: $names_hash{$parent} ($nodes_hash{$parent})\n";
      }
    }
    return $result;
}

sub tax_walk_to_get_rank {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
  #print "$parent, $parent_rank\n";
  my ($phylum,$kingdom,$superkingdom) = ("undef","undef","undef");
  my $result;

  while (1) {
    if ($parent_rank eq "phylum") {
      $phylum = $names_hash{$parent};
      #print "Found phylum: $phylum\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "kingdom") {
      $kingdom = $names_hash{$parent};
      #print "Found phylum: $kingdom\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "superkingdom") {
      $superkingdom = $names_hash{$parent};
      #print "Found phylum: $superkingdom\n";
      last;
    } elsif ($parent == 1) {
      last;
    } else {
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
    }
  }
  $result = join (";",$superkingdom,$kingdom,$phylum);
  return $result;
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

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
