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
  -f|--fofn       [FILE]   : fofn of treefiles, one per line
  -o|--out_suffix [FILE]   : suffix to be added to modified treefile ['.tax.treefile']
  -p|--taxdb      [STRING] : path to nodes.dmp and names.dmp tax files
  -s|--speclist   [FILE]   : UniProt species identification codes file (https://www.uniprot.org/docs/speclist)
  -t|--taxlist    [FILE]   : UniProt taxid file, formatted 'uniprotid TAB taxid' (optional)
  -d|--taxdepth   [STRING] : depth of taxonomy returned, can be 'phylum' or 'species' ['phylum']
  -n|--taxnumber           : include taxid integer in string? [no]
  -h|--help                : prints this help message
\n";

my ($infile,$fofn,$path,$tax_list,$spec_list,$tax_number,$help);
my $out_suffix = "tax.treefile";
my $tax_depth = "phylum";

GetOptions (
  'i|infile:s'  => \$infile,
  'f|fofn:s'  => \$fofn,
  'o|out_suffix:s' => \$out_suffix,
  'p|taxdb=s'   => \$path,
  's|speclist=s' => \$spec_list,
  't|taxlist:s' => \$tax_list,
  'd|depth:i'   => \$tax_depth,
  'n|taxnumber' => \$tax_number,
  'h|help'      => \$help,
);

die $usage if $help;
die $usage unless (($infile || $fofn) && $spec_list && $path);

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

my @treefiles;
if ( $infile ) {
  push (@treefiles, $infile);
  print STDERR "[INFO] Working with one treefile: '$infile'\n";
} elsif ( $fofn ) {
  open (my $FOFN, $fofn) or die $!;
  while (my $line = <$FOFN>) {
    chomp $line;
    push (@treefiles, $line);
  }
  close $FOFN;
  print STDERR "[INFO] Number of treefiles in '$fofn': ".scalar(@treefiles)."\n";
}

foreach my $current_file (@treefiles) {

  my %tax_hash;
  my @uniprot_strings;

  open (my $TREEFILE_READ, $current_file) or die $!;
  while (<$TREEFILE_READ>) {
    ## this regex *should* capture most UniProt accessions and gene IDs, plus the (_\d+\-\d+){0,1} suffix added by hmmalign/IQTREE
    while ($_ =~ m/(\w{1,5}_\w{1,5}(?<=[A-Z])(_\d+\-\d+){0,1}|[OPQ][0-9][A-Z0-9]{3}[0-9]_[A-Z0-9]{1,5}(_\d+\-\d+){0,1}|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}_[A-Z0-9]{1,5}(_\d+\-\d+){0,1})/g) {
      push (@uniprot_strings, $1);
    }
    ## explanation
    ## m/(\w{1,5}_\w{1,5}(?<=[A-Z]) ## match gene ID followed by species ID, requires at least one letter e.g. CASP2_RAT
    ## (_\d+\-\d+){0,1} ## maybe match hmmalign stuff
    ## |[OPQ][0-9][A-Z0-9]{3}[0-9]_[A-Z0-9]{1,5}(_\d+\-\d+){0,1}|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}_[A-Z0-9]{1,5}(_\d+\-\d+){0,1})/g ## or match UniProt ID
  }
  close $TREEFILE_READ;
  print STDERR "[INFO] Number of UniProt IDs scooped from treefile '$current_file': ".commify(scalar(@uniprot_strings))."\n";

  foreach my $orig_string (@uniprot_strings) {
    ## default is the original ID
    my $taxid;
    my @a = split ("_", $orig_string);
    my $replace_string = join ("_", $a[0], $a[1]); ## default is to include the UniProt species code
    ## try to grep from speclist using species ID first for speed
    my $match = `grep -m1 -wF $a[1] $spec_list`;
    if ( $match ) {
      my @b = split (m/\s+/, $match);
      $taxid = $b[2];
      $taxid =~ s/://;
    } else {
      ## else parse UniProt accession for taxid
      $match = `grep -m1 -wF $a[0] $tax_list`;
      my @b = split (m/\s+/, $match);
      $taxid = $b[1];
    }
    ## pull taxonomy string from taxdb
    if ( $taxid ) {
      if (($taxid =~ m/\d+/) && (check_taxid_has_parent($taxid) == 0)) {
        if ($tax_depth eq "species") {
          print STDERR " --> " . join (" ", join("_",$a[0],$a[1]), $taxid, tax_walk_to_get_rank_to_species($taxid)) . "\n";
          if ( $tax_number ) {
            $replace_string = join ("_", $a[0], $a[1], tax_walk_to_get_rank_to_species($taxid), $taxid);
          } else {
            $replace_string = join ("_", $a[0], $a[1], tax_walk_to_get_rank_to_species($taxid));
          }
        } else {
          print STDERR " --> " . join (" ", join("_",$a[0],$a[1]), $taxid, tax_walk_to_get_rank_to_phylum($taxid)) . "\n";
          if ( $tax_number ) {
            $replace_string = join ("_", $a[0], $a[1], tax_walk_to_get_rank_to_phylum($taxid), $taxid);
          } else {
            $replace_string = join ("_", $a[0], $a[1], tax_walk_to_get_rank_to_phylum($taxid));
          }
        }
      } else {
        print STDERR " --> " . join (" ", join("_",$a[0],$a[1]), "Invalid taxid") . "\n";
      }
    } else {
      print STDERR " --> " . join (" ", join("_",$a[0],$a[1]), "No taxid found") . "\n";
    }
    ## replacement string hash
    $tax_hash{$orig_string} = $replace_string;
  }

  ## set up regex
  my $regex = join ("|", keys %tax_hash);
  $regex = qr/$regex/;

  ## open treefile again and make the substitution:
  print STDERR "[INFO] Printing new tree to '$current_file.$out_suffix'...\n";

  open (my $TREEFILE_WRITE, $current_file) or die $!;
  open (my $OUT, ">$current_file.$out_suffix") or die $!;
  while (my $tree = <$TREEFILE_WRITE>) {
    $tree =~ s/($regex)/$tax_hash{$1}/g;
    print $OUT $tree;
  }
  close $TREEFILE_WRITE;
  close $OUT;
}

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
  $result =~ s/_$//; ## remove trailing underscore
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
  $result =~ s/_$//; ## remove trailing underscore
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
