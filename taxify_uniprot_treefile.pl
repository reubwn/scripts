#!/usr/bin/env perl

## author: reubwn April 2021

use strict;
use warnings;

use Getopt::Long;

my $usage = "
SYNOPSIS:
  Annotates a treefile containing UniRef sequence IDs with taxonomy information.
  Specifically, it replaces the mnemonic UniProt species code with a tax string (and optional taxid) to either species or phylum (default) level.

OPTIONS:
  -i|--infile     [FILE]   : input treefile
  -o|--out_suffix [FILE]   : suffix to be added to modified treefile ['.tax.treefile']
  -p|--taxdb      [STRING] : path to nodes.dmp and names.dmp tax files
  -t|--taxids     [FILE]   : UniProt taxid file, formatted 'uniprotid TAB taxid'
  -d|--taxdepth   [STRING] : depth of taxonomy returned, can be 'phylum' or 'species' ['phylum']
  -n|--taxnumber           : include taxid integer in string? [no]
  -h|--help                : prints this help message
\n";

my ($infile,$path,$tax_list,$tax_number,$help);
my $out_suffix = "tax.treefile";
my $tax_depth = "phylum";

GetOptions (
  'i|infile=s'  => \$infile,
  'o|out_suffix:s' => \$out_suffix,
  'p|taxdb=s'   => \$path,
  't|taxids=s'  => \$tax_list,
  'd|depth:i'   => \$tax_depth,
  'n|taxnumber' => \$tax_number,
  'h|help'      => \$help,
);

die $usage if $help;
die $usage unless ($infile && $tax_list && $path);

## parse nodes and names:
my (%nodes_hash, %names_hash, %rank_hash);

print STDERR "[INFO] Building taxonomy databases from tax files in '$path'...\n";
open (my $NODES, "$path/nodes.dmp") or die $!;
while (<$NODES>) {
  chomp;
  next if /\#/;
  my @F = map { s/^\s+|\s+$//gr } split (m/\|/, $_); ## split nodes.dmp file on \s+|\s+ regex
  $nodes_hash{$F[0]} = $F[1]; ## key= child taxid; value= parent taxid
  $rank_hash{$F[0]} = $F[2]; ## key= taxid; value= rank
}
close $NODES;
open (my $NAMES, "$path/names.dmp") or die $!;
while (<$NAMES>) {
  chomp;
  next if /\#/;
  my @F = map { s/^\s+|\s+$//gr } split (m/\|/, $_);
  $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
}
close $NAMES;
if (-e "$path/merged.dmp") {
  open (my $MERGED, "$path/merged.dmp") or die $!;
  while (<$MERGED>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (m/\|/, $_);
    $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
    ## this will behave as if old taxid is a child of the new one, which is OK I guess
  }
}
print STDERR "[INFO] Nodes parsed: ".commify(scalar(keys %nodes_hash))."\n";

## parse treefile:
print STDERR "[INFO] Getting UniProt IDs from '$infile'...\n";

my %tax_hash;
open (my $TREEFILE_READ, $infile) or die $!;
while (<$TREEFILE_READ>) {
  ## regex to capture UniProt ID string
  my @uniprot_strings = ($_ =~ m/([OPQ][0-9][A-Z0-9]{3}[0-9]_[A-Z0-9]{1,5}_\d+\-\d+|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]{1,2}_[A-Z0-9]{1,5}_\d+\-\d+)/g);
  print STDERR "[INFO] Number of UniProt IDs in treefile: " . commify(scalar(@uniprot_strings)) . "\n";
  foreach my $orig_string (@uniprot_strings) {
    my @a = split ("_", $orig_string);
    my $replace_string = join ("_", $a[0], $a[1]); ## default is to include the UniProt species code
    ## grep UniProt ID from UniProt taxid file
    my $match = `grep -m1 -wF $a[0] $tax_list`;
    my @b = split (m/\s+/, $match);
    if (($b[1] =~ m/\d+/) && (check_taxid_has_parent($b[1]) == 0)) {
      if ($tax_depth eq "species") {
        print STDERR " --> " . join (" ", join("_",$a[0],$a[1]), $b[1], tax_walk_to_get_rank_to_species($b[1])) . "\n";
        if ( $tax_number ) {
          $replace_string = join ("_", $a[0], tax_walk_to_get_rank_to_species($b[1]), $b[1]);
        } else {
          $replace_string = join ("_", $a[0], tax_walk_to_get_rank_to_species($b[1]));
        }
      } else {
        print STDERR " --> " . join (" ", join("_",$a[0],$a[1]), $b[1], tax_walk_to_get_rank_to_phylum($b[1])) . "\n";
        if ( $tax_number ) {
          $replace_string = join ("_", $a[0], tax_walk_to_get_rank_to_phylum($b[1]), $b[1]);
        } else {
          $replace_string = join ("_", $a[0], tax_walk_to_get_rank_to_phylum($b[1]));
        }
      }
    } else {
      print STDERR " --> " . join (" ", join("_",$a[0],$a[1]), "Invalid TaxID") . "\n";
    }
    ## replacement string hash
    $tax_hash{$orig_string} = $replace_string;
  }
}
close $TREEFILE_READ;

## set up regex
my $regex = join ("|", keys %tax_hash);
$regex = qr/$regex/;

## open treefile again and make the substitution:
print STDERR "[INFO] Printing new tree to '$infile.$out_suffix'...\n";

open (my $TREEFILE_WRITE, $infile) or die $!;
open (my $OUT, ">$infile.$out_suffix") or die $!;
while (my $tree = <$TREEFILE_WRITE>) {
  $tree =~ s/($regex)/$tax_hash{$1}/g;
  print $OUT $tree;
}
close $TREEFILE_WRITE;
close $OUT;

print STDERR "[INFO] Done " . `date`;

###### SUBS

sub check_taxid_has_parent {
  my $taxid = $_[0];
  my $result = 0;
  unless ($nodes_hash{$taxid}) {
    $result = 1;
  }
  return $result; ## 0 = taxid exists; 1 = taxid does not exist
}

sub tax_walk_to_get_rank_to_phylum {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
  my ($superkingdom,$kingdom,$phylum,) = ("undef","undef","undef");

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
  my $result = join ("_",$superkingdom,$kingdom,$phylum);
  $result =~ s/\s+/\_/g; ## replace spaces with underscores
  $result =~ s/undef//g; ## remove undef
  $result =~ s/__/_/g; ## replace any __ with _
  return $result;
}

sub tax_walk_to_get_rank_to_species {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
  my ($superkingdom,$kingdom,$phylum,$class,$order,$family,$genus,$species) = ("undef","undef","undef","undef","undef","undef","undef","undef");

  while (1) {
    if ($parent_rank eq "species") {
      $species = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "genus") {
      $genus = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "family") {
      $family = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "order") {
      $order = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "class") {
      $class = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "phylum") {
      $phylum = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "kingdom") {
      $kingdom = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "superkingdom") {
      $superkingdom = $names_hash{$parent};
      last;
    } elsif ($parent == 1) {
      last;
    } else {
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
    }
  }
  my $result = join ("_",$superkingdom,$kingdom,$phylum,$class,$order,$family,$genus,$species);
  $result =~ s/\s+/_/g; ## replace spaces with underscores
  $result =~ s/undef//g; ## remove undef
  $result =~ s/__/_/g; ## replace any __ with _
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
