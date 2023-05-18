#!/usr/bin/env perl

## reubwn May 23

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Term::ANSIColor;
use File::Path qw( rmtree );
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Convert OrthoFinder orthogroups file (N0.tsv or Orthogroups.txt) to phyletic presence/absence matrix.
  Globs '*.fasta', '*.faa', '*.fa' or '*.aa' from --path.

OPTIONS [*required]
  -i|--groups  *[FILE] : Orthogroups file (N0.tsv or Orthogroups.txt)
  -p|--path    *[PATH] : path to protein files
  -l|--legacy          : file is Orthogroups.txt format (default is N0.tsv format)
  -u|--unassigned      : add counts of unassigned genes (lineage-specific) to matrix (default: don't add)
  -o|--out      [STR]  : outfiles prefix ('phyletic')
  -h|--help            : print this message
\n";

my ($orthogroups_file, $proteins_path, $legacy, $add_unassigned, $help);
my $outprefix = "out";

GetOptions (
  'i|groups=s' => \$orthogroups_file,
  'p|path=s'   => \$proteins_path,
  'l|legacy'   => \$legacy,
  'u|unassigned' => \$add_unassigned,
  'o|out:s'    => \$outprefix,
  'h|help'     => \$help
);

die $usage if ( $help );
die $usage unless ( $orthogroups_file && $proteins_path );

## get sequences
my %seq_hash;
my @fastas = glob("$proteins_path/*.fasta $proteins_path/*.faa $proteins_path/*.fa $proteins_path/*.aa");
if (scalar(@fastas) == 0) {
  print STDERR "[ERROR] Nothing found in '$proteins_path'\n";
  print STDERR "[ERROR] Please make sure there are protein fastas in '$proteins_path' with extension fasta|faa|fa|aa\n";
  die 1;
} else {
  print STDERR "[INFO] Found ".scalar(@fastas)." files in '$proteins_path'\n";
}

print STDERR "[INFO] Reading sequences from:\n";
my %species_hash;
foreach (@fastas){
  print STDERR colored($_, 'white on_blue') . "\n";
  my $in = Bio::SeqIO->new ( -file => $_, -format => "fasta" );
  while ( my $seq_obj = $in->next_seq() ){
    my $seq_id = $seq_obj->display_id();
    my $species_id = ( split(m/\|/, $seq_id) )[0];
    $seq_hash{$seq_id} = ();
    $species_hash{$species_id} = ();
  }
}
print STDERR "[INFO] Read in ".commify(scalar(keys %seq_hash))." sequences from ".commify(scalar(@fastas))." files\n";

## open Orthogroups
my %orthogroups;
my %seen_in_orthogroups;
my %unassigned_seqs = %seq_hash;

# matrix looks like
# >species1
# 00010101100111101101111
# where cols are OGs and membership is collapsed to 1/0

open (my $GROUPS, $orthogroups_file) or die $!;
while (my $line = <$GROUPS>) {
  chomp $line;
  my @seqs_array;

  if ( $legacy ) { ## Orthogroups.txt format
    # print STDERR "[INFO] Input file is 'Orthogroups.txt' format\n";
    my @a = split(/\:\s+/, $line);
    my @b = split(/\s+/, $a[1]);
    @seqs_array = @b;
    print STDERR "\r[INFO] Working on OG \#$.: $a[0]"; $|=1;

  } else { ## N0.tsv format
    # print STDERR "[INFO] Input file is 'N0.tsv' (HOGs) format\n";
    next if $. == 1; ## skip first line
    my @a = split(/\s+/, $line);
    my @b = @a[3..$#a]; ## seqs begin in column 4
    my @c = map{s/,//g; $_} @b;
    @seqs_array = @c;
    print STDERR "\r[INFO] Working on HOG \#$.: $a[0]"; $|=1;
  }

  my %membership_per_OG_hash;
  ## collapse OG membership to 1/0
  foreach (@seqs_array) {
    my $species_id = ( split(/\|/, $_) )[0];
    $membership_per_OG_hash{$species_id}++;

    ## delete from %unassigned_seqs; only genes not found in OGs remain in this hash
    delete $unassigned_seqs{$_};
  }

  ## cycle thru species hash and push 1/0 depending on membership
  foreach (nsort keys %species_hash) {
    if ($membership_per_OG_hash{$_}) {
      push ( @{$species_hash{$_}}, 1 );
    } else {
      push ( @{$species_hash{$_}}, 0 );
    }
  }
}
close $GROUPS;

my %unassigned_counts;
foreach (keys %unassigned_seqs) {
  my @a = split (/\|/, $_);
  $unassigned_counts{$a[0]}++; ## should be the number of seqs that weren't found in any HOG
}

print STDERR "\n[INFO] Unassigned gene counts:\n";
foreach (nsort keys %unassigned_counts) {
  print STDERR "$_ : $unassigned_counts{$_}\n";
}

## add lineage specific unassigned genes to the end of matrix
if ($add_unassigned) {
  foreach my $k1 (nsort keys %species_hash) {
    foreach my $k2 (nsort keys %unassigned_counts) {
      if ($k1 eq $k2) {
        ## adds 1 to correct species to length == number of unassigned genes
        for (1..$unassigned_counts{$k2}) {
          push ( @{$species_hash{$k1}}, 1 );
        }
      } else {
        ## else adds 0 to any other species
        for (1..$unassigned_counts{$k2}) {
          push ( @{$species_hash{$k1}}, 0 );
        }
      }
    }
  }
}

## open outfile and print as fasta phyletic matrix
open (my $OUT, ">$outprefix".".phyletic_matrix.txt");
foreach (nsort keys %species_hash) {
  print $OUT ">$_\n";
  print $OUT join ('', @{$species_hash{$_}}) . "\n";
}
close $OUT;

print STDERR "[INFO] Printing matrix to file '$outprefix.phyletic_matrix.txt'\n";
print STDERR "[INFO] Finished on ".`date`."\n";

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
