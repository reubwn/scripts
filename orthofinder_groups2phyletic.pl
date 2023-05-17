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
  Convert OrthoFinder Orthogroups.txt file to phyletic presence/absence matrix.
  Globs '*.fasta', '*.faa', '*.fa' or '*.aa' from --path.

OPTIONS [*required]
  -i|--groups  *[FILE] : Orthogroups.txt file
  -p|--path    *[PATH] : path to protein files
  -o|--out      [STR]  : outfiles prefix ('phyletic')
  -h|--help            : print this message
\n";

my ($orthogroups_file, $proteins_path, $help);
my $outprefix = "out";

GetOptions (
  'i|groups=s' => \$orthogroups_file,
  'p|path=s'   => \$proteins_path,
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
print STDERR "[INFO] Read in ".commify(scalar(keys %seq_hash))." sequences from ".commify(scalar(@fastas))." files\n\n";
# foreach (keys %species_hash) {
#   print "$_\n";
# }

## open Orthogroups
my %orthogroups;
my %seen_in_orthogroups;

# matrix looks like
# >species1
# 00010101100111101101111
# where cols are OGs and membership is collapsed to 1/0

open (my $GROUPS, $orthogroups_file) or die $!;

while (my $line = <$GROUPS>) {
  chomp $line;
  my @a = split(/\:\s+/, $line);
  my @b = split(/\s+/, $a[1]);
  print STDERR "\r[INFO] Working on OG \#$.: $a[0]"; $|=1;

  my %membership_per_OG_hash;
  ## collapse OG membership to 1/0
  foreach (@b) {
    my $species_id = ( split(/\|/, $_) )[0];
    $membership_per_OG_hash{$species_id}++;
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

## open outfile and print as fasta phyletic matrix
open (my $OUT, ">$outprefix".".phyletic_matrix.txt");
foreach (nsort keys %species_hash) {
  print $OUT ">$_\n";
  print $OUT join ('', @{$species_hash{$_}}) . "\n";
}
close $OUT;

print STDERR "\n[INFO] Finished on ".`date`."\n";




######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
