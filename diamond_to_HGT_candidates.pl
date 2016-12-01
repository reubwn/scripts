#!/usr/bin/env perl

## author: reubwn Nov 2016

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use Data::Dumper qw(Dumper);
use List::Util qw(reduce sum);

my $usage = "
SYNOPSIS:
  Goal is to take a taxified Diamond or BLAST file, and for each hit recurse up the
  tax tree until that hit can be categorised into \"ingroup\" versus \"outgroup\"
  (e.g., \"metazoan\" vs \"non-metazoan\" etc.).
    (a) Get Query Taxid: For each query, calculate the sum of bitscores per taxid;
        the taxid with the highest bitscoresum is the \"winner\" and assumed taxonomy of the query
    (b) Determine Category: Assess the placement of the winning taxid in the NCBI
        taxonomy tree: if the taxid falls within the --taxid_threshold (default = Metazoa)
        then the query is determined to be INGROUP; if the taxid does not fall within the
        taxid_threshold then the query is OUTGROUP
    (c) Get Support: Assess support for the winning query taxid from secondary hits;
        winning taxid is well-supported if the rest of the hits agree with the INGROUP/OUTGROUP
        categorisation above --support_threshold (default = 90%)
    (d) Print: INGROUP/OUTGROUP decisions and support values are printed to file for all
        queries, regardless of support

    Eg: The bestsum bitscore for a given query protein is to a Pseudomonad species; thus taxid is
    ___ and category is OUTGROUP as taxid ___ does not fall within Metazoa. The vast majority (>90\%)
    of other hits are also to bacteria, and thus support the OUTGROUP categorisation. This gene
    would be designated \"candidate HGT\" and would merit further investigation.

OUTPUTS:
  A \"\*.HGT_decisions.\" file with the headers:

OPTIONS:
  -i|--in                [FILE]   : tab formatted Diamond output file [required]
  -p|--path              [STRING] : path to dir/ containing nodes.dmp, names.dmp and merged.dmp
  -o|--nodes             [FILE]   : path to nodes.dmp [required unless --nodesDB]
  -a|--names             [FILE]   : path to names.dmp [required unless --nodesDB]
  -m|--merged            [FILE]   : path to merged.dmp
  -n|--nodesDB           [FILE]   : nodesDB.txt file from blobtools [required unless --nodes && --names]
  -t|--taxid_threshold   [INT]    : NCBI taxid to recurse up to; i.e., threshold taxid to define \"ingroup\" [default = 33208 (Metazoa)]
  -s|--support_threshold [FLOAT]  : support required to be counted as \"well-supported\" [default = 90\%; note all are printed regardless of support]
  -b|--bitscore_column   [INT]    : define bitscore column for --in (first column = 1) [default: 12]
  -c|--taxid_column      [INT]    : define taxid column for --in (first column = 1) [default: 13]
  -d|--delimiter         [STRING] : define delimiter to split --in (specify \"diamond\" for Diamond files (\"\\s+\") or \"blast\" for BLAST files (\"\\t\")) [default: diamond]
  -e|--prefix            [FILE]   : filename prefix for outfile [default = INFILE]
  -H|--header                     : don't print header [default: do print it]
  -v|--verbose                    : say more things [default: be quiet]
  -h|--help                       : prints this help message

EXAMPLES:

\n";

my ($in,$nodesfile,$path,$namesfile,$mergedfile,$nodesDBfile,$prefix,$outfile,$hgtcandidatesfile,$warningsfile,$header,$verbose,$debug,$help);
my $taxid_threshold = 33208;
my $support_threshold = 90;
my $taxid_column = 13;
my $bitscore_column = 12;
my $delimiter = "diamond";

GetOptions (
  'in|i=s'              => \$in,
  'path|p:s'            => \$path,
  'nodes|o:s'           => \$nodesfile,
  'names|a:s'           => \$namesfile,
  'merged|m:s'          => \$mergedfile,
  'nodesDB|n:s'         => \$nodesDBfile,
  'taxid_threshold|t:i' => \$taxid_threshold,
  'support_threshold|s:f' => \$support_threshold,
  'taxid_column|c:i'    => \$taxid_column,
  'bitscore_column|b:i' => \$bitscore_column,
  'delimiter|d:s'       => \$delimiter,
  'prefix|e:s'          => \$prefix,
  'header|H'            => \$header,
  'verbose|v'           => \$verbose,
  'debug'               => \$debug,
  'help|h'              => \$help,
);

die $usage if $help;
die $usage unless ($in);
#die $usage unless (($nodesfile && $namesfile) || $nodesDBfile || $path);

## define delimiter:
if ($delimiter eq "diamond") {
  $delimiter = qr/\s+/;
} elsif ($delimiter eq "blast") {
  $delimiter = qr/\t/;
} else {
  die "[ERROR] Unknown delimiter, please choose \"diamond\" or \"blast\"\n";
}

############################################## PARSE NODES

## parse nodes and names:
my (%nodes_hash, %names_hash, %rank_hash);
if ($path) {
  print STDOUT "[INFO] Building taxonomy databases from tax files in \"$path\"...";
  open(my $NODES, "$path/nodes.dmp") or die $!;
  while (<$NODES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_); ## split nodes.dmp file on \s+|\s+ regex
    $nodes_hash{$F[0]} = $F[1]; ## key= child taxid; value= parent taxid
    $rank_hash{$F[0]} = $F[2]; ## key= taxid; value= rank
  }
  close $NODES;
  open (my $NAMES, "$path/names.dmp") or die $!;
  while (<$NAMES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
  }
  close $NAMES;
  if (-e "$path/merged.dmp") {
    open (my $MERGED, "$path/merged.dmp") or die $!;
    while (<$MERGED>) {
      chomp;
      next if /\#/;
      my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
      $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
      ## this will behave as if old taxid is a child of the new one, which is OK I guess
    }
  }
} elsif ($nodesfile && $namesfile) {
  print STDOUT "[INFO] Building taxonomy databases from \"$nodesfile\" and \"$namesfile\"...";
  open(my $NODES, $nodesfile) or die $!;
  while (<$NODES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_); ## split nodes.dmp file on \s+|\s+ regex
    $nodes_hash{$F[0]} = $F[1]; ## key= child taxid; value= parent taxid
    $rank_hash{$F[0]} = $F[2]; ## key= taxid; value= rank
  }
  close $NODES;
  open (my $NAMES, $namesfile) or die $!;
  while (<$NAMES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
  }
  close $NAMES;
  if ($mergedfile) {
    open (my $MERGED, $mergedfile) or die $!;
    while (<$MERGED>) {
      chomp;
      next if /\#/;
      my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
      $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
      ## this will behave as if old taxid is a child of the new one, which is OK I guess
    }
  }
} elsif ($nodesDBfile) {
  print STDOUT "[INFO] Building taxonomy databases from \"$nodesDBfile\"..."; $| = 1;
  open(my $NODES, $nodesDBfile) or die $!;
  while (<$NODES>) {
    chomp;
    next if /\#/;
    my @F = split (/\t/, $_);
    $nodes_hash{$F[0]} = $F[3]; ## key= child taxid; value= parent taxid
    $names_hash{$F[0]} = $F[2]; ## key= taxid; value= species name
    $rank_hash{$F[0]} = $F[1]; ## key= taxid; value= rank
  }
  close $NODES;
}
## print some info to STDOUT:
print STDOUT " done\n";
print STDOUT "[INFO] Nodes parsed: ".scalar(keys %nodes_hash)."\n";
print STDOUT "[INFO] Threshold taxid set to \"$taxid_threshold\" ($names_hash{$taxid_threshold})\n";
print STDOUT "[INFO] INGROUP set to \"$names_hash{$taxid_threshold}\"; OUTGROUP is therefore \"non-$names_hash{$taxid_threshold}\"\n";

############################################ OUTFILES

## define outfiles:
if ($prefix) {
  $outfile = "$prefix.HGT_decisions.$names_hash{$taxid_threshold}.txt";
  $hgtcandidatesfile = "$prefix.HGT_candidates.$names_hash{$taxid_threshold}.$support_threshold.txt";
  $warningsfile = "$prefix.warnings.txt";
} else {
  $outfile = "$in.HGT_decisions.$names_hash{$taxid_threshold}.txt";
  $hgtcandidatesfile = "$in.HGT_candidates.$names_hash{$taxid_threshold}.$support_threshold.txt";
  $warningsfile = "$in.warnings.txt";
}

## open outfiles:
open (my $OUT, ">",$outfile) or die $!;
open (my $HGT, ">",$hgtcandidatesfile) or die $!;
open (my $WARN, ">",$warningsfile) or die $!;
print $OUT "\#query\ttaxid\tbestsumBitscore\tsuperkingdom;kingdom;phylum\tingroupTaxname\tdecision\tsupport\n" unless $header;
print $HGT "\#query\ttaxid\tbestsumBitscore\tsuperkingdom;kingdom;phylum\tingroupTaxname\tdecision\tsupport\n" unless $header;

############################################## PARSE DIAMOND

## parse Diamond file:
print STDOUT "[INFO] Parsing Diamond file \"$in\"...";
my (%sum_bitscores_per_query_hash);
my $skipped_entries = 0;
open (my $DIAMOND, $in) or die $!;
while (<$DIAMOND>) {
  chomp;
  next if /^\#/;
  my @F = split ($delimiter, $_);
  if ($F[($taxid_column-1)] !~ m/\d+/) {
    print $WARN "[WARN] The taxid ".$F[($taxid_column-1)]." for query $F[0] on line $. of \"$in\" does not look like a valid NCBI taxid... Skipping this entry\n";
    $skipped_entries++;
    next;
  } elsif (check_taxid_has_parent($F[($taxid_column-1)]) == 1) {
    print $WARN "[WARN] The taxid ".$F[($taxid_column-1)]." for query $F[0] on line $. of \"$in\" does not appear to have a valid parent tax node... Skipping this entry\n";
    $skipped_entries++;
    next;
  }

  ## push all bitscores for every taxid into an array within a hash within a hash:
  push @{ $sum_bitscores_per_query_hash{$F[0]}{$F[($taxid_column-1)]} }, $F[($bitscore_column-1)]; ## key= query; value = hash{ key= taxid; value= array[ bitscores ]}
}
close $DIAMOND;
print STDOUT " done\n";
print STDOUT "[WARN] There were $skipped_entries invalid taxids in \"$in\"; these entries were omitted from the analysis.\n" if $skipped_entries > 0;

############################################ DEBUG

print Dumper \%sum_bitscores_per_query_hash if $debug;

############################################ MAIN CODE

## get winning bitscore and taxid; calculate congruence among all taxids for all hits per query:
my ($processed,$ingroup_supported,$outgroup_supported) = (0,0,0);
print STDOUT "[INFO] Calculating bestsum bitscore and hit support...\n";
print STDOUT "\n" if $verbose;
foreach my $query (nsort keys %sum_bitscores_per_query_hash) {
  my %bitscore_hash = %{ $sum_bitscores_per_query_hash{$query} }; ## key= taxid; value= \@array of all bitscores for that taxid

  ## calculate bitscoresums per taxid; get taxid of highest bitscoresum; get support for winning taxid from other hits:
  my (%bitscoresum_hash, %count_categories, %support_categories);
  foreach my $taxid (nsort keys %bitscore_hash) {
    my $bitscoresum = sum( @{ $bitscore_hash{$taxid} } );
    $bitscoresum_hash{$taxid} = $bitscoresum; ## key= taxid; value= bitscoresum
    $count_categories{tax_walk($taxid)}++; ## count categories; if each hit's taxid falls within/outwith the $taxid_threshold
  }
  print "$query:\n" if $debug; ## debug
  print Dumper \%bitscoresum_hash if $debug; ## debug

  foreach my $cat (keys %count_categories) {
    $support_categories{$cat} = percentage($count_categories{$cat}, scalar(keys %bitscore_hash)); ## calculate proportion of support for the category of the winner
  }

  ## get taxid with highest bitscore:
  my $taxid_with_highest_bitscore = List::Util::reduce { $bitscoresum_hash{$b} > $bitscoresum_hash{$a} ? $b : $a } keys %bitscoresum_hash; ## winning taxid
  my $taxid_with_highest_bitscore_category = tax_walk($taxid_with_highest_bitscore); ## category of winning taxid ("ingroup", "outgroup" or "unassigned")
  my $taxid_with_highest_bitscore_category_support = $support_categories{$taxid_with_highest_bitscore_category}; ## % support from other hits
  print "[INFO] [$query] Taxid with highest bitscore: $taxid_with_highest_bitscore (bitscore = $bitscoresum_hash{$taxid_with_highest_bitscore}; taxonomy = ".tax_walk_to_get_rank_to_phylum($taxid_with_highest_bitscore).")\n" if $verbose;
  print "[INFO] [$query] Decision of bestsum bitscore: $taxid_with_highest_bitscore_category (support = $taxid_with_highest_bitscore_category_support)\n" if $verbose;

  ## print to $out
  print $OUT join "\t", $query, $taxid_with_highest_bitscore, $bitscoresum_hash{$taxid_with_highest_bitscore}, tax_walk_to_get_rank_to_phylum($taxid_with_highest_bitscore), "ingroup=".$names_hash{$taxid_threshold}, $taxid_with_highest_bitscore_category, $taxid_with_highest_bitscore_category_support, "\n";

  ## calculate number of well-supported genes:
  if ( ($taxid_with_highest_bitscore_category eq "ingroup") && ($taxid_with_highest_bitscore_category_support > $support_threshold) ) {
    $ingroup_supported++;
  } elsif ( ($taxid_with_highest_bitscore_category eq "outgroup") && ($taxid_with_highest_bitscore_category_support > $support_threshold) ) {
    print $HGT join "\t", $query, $taxid_with_highest_bitscore, $bitscoresum_hash{$taxid_with_highest_bitscore}, tax_walk_to_get_rank_to_species($taxid_with_highest_bitscore), "ingroup=".$names_hash{$taxid_threshold}, $taxid_with_highest_bitscore_category, $taxid_with_highest_bitscore_category_support, "\n";
    $outgroup_supported++;
  }

  ## progress
  $processed++;
  if ($processed % 1000 == 0){
    print "\r[INFO] Processed ".commify($processed)." queries...";
    $| = 1;
  }
}
close $OUT;
close $HGT;
close $WARN;
print STDOUT "\r[INFO] Processed ".commify($processed)." queries\n";
print STDOUT "[INFO] Number of queries in category \"$names_hash{$taxid_threshold}\" with support > $support_threshold\%: ".commify($ingroup_supported)."\n";
print STDOUT "[INFO] Number of queries in category \"non-$names_hash{$taxid_threshold}\" with support > $support_threshold\%: ".commify($outgroup_supported)."\n";
print STDOUT "[INFO] Finished\n\n";

############################################ SUBS

sub check_taxid_has_parent {
  my $taxid = $_[0];
  my $result = 0;
  unless ($nodes_hash{$taxid}) {
    $result = 1;
  }
  return $result; ## 0 = taxid exists; 1 = taxid does not exist
}

sub tax_walk {
    my $taxid = $_[0];
    my $walk_to;
    if (exists $_[1]) {
      $walk_to = $_[1];
    } else {
      $walk_to = $taxid_threshold; ## default is metazoa
    }

    ## first parent:
    my $parent = $nodes_hash{$taxid};
    my $result;

    ## return "unassigned" if hit has no valid taxid
    if ($parent !~ m/\d+/) {
      $result = "unassigned";
      return $result;
    }

    ## recurse the tree:
    while (1) {
      if ($parent == $walk_to) {
        $result = "ingroup"; ## is ingroup
        last;
      } elsif ($parent == 1) {
        $result = "outgroup"; ## root; i.e., the whole tree has been recursed without finding $threshold, therefore $taxid must reside in another part of the tax tree
        last;
      } elsif ($parent == 32644) {
        $result = "unassigned"; ## taxid for "unidentified"
        last;
      } else { ## walk up the tree!
        $parent = $nodes_hash{$parent};
      }
    }
    return $result;
}

sub tax_walk_to_get_rank_to_phylum {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
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

sub tax_walk_to_get_rank_to_species {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
  my ($species,$genus,$family,$order,$class,$phylum,$kingdom,$superkingdom) = ("undef","undef","undef","undef","undef","undef","undef","undef");
  my $result;

  while (1) {
    if ($parent_rank eq "species") {
      $phylum = $names_hash{$parent};
      #print "Found phylum: $phylum\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "genus") {
      $kingdom = $names_hash{$parent};
      #print "Found phylum: $kingdom\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "family") {
      $kingdom = $names_hash{$parent};
      #print "Found phylum: $kingdom\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "order") {
      $kingdom = $names_hash{$parent};
      #print "Found phylum: $kingdom\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "class") {
      $kingdom = $names_hash{$parent};
      #print "Found phylum: $kingdom\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "phylum") {
      $kingdom = $names_hash{$parent};
      #print "Found phylum: $kingdom\n";
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
  $result = join (";",$superkingdom,$kingdom,$phylum,$class,$order,$family,$genus,$species);
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
