#!/usr/bin/env perl

## reubwn May 23

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Convert OrthoFinder N0.tsv file to legacy Orthogroups.txt format (OrthoMCL-style) for downstream tools.
  Finds unassigned genes directly from N0.tsv and protein files and adds them as 'singleton OGs' at the end of the file.
  Globs '*.fasta', '*.faa', '*.fa' or '*.aa' from --path.

OPTIONS [*required]
  -i|--groups  *[FILE] : OrthoFinder N0.tsv file
  -p|--path    *[PATH] : path to protein files
  -U|--no_unassigned   : DON'T add counts of unassigned genes to matrix (default: add them)
  -o|--out      [STR]  : outfiles prefix ('phyletic')
  -h|--help            : print this message
\n";

my ($orthogroups_file, $proteins_path, $no_unassigned, $help);
my $outprefix = "out";

GetOptions (
  'i|groups=s' => \$orthogroups_file,
  'p|path=s'   => \$proteins_path,
  'U|no_unassigned' => \$no_unassigned,
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
  print STDERR "\t$_\n";
  my $in = Bio::SeqIO->new ( -file => $_, -format => "fasta" );
  while ( my $seq_obj = $in->next_seq() ){
    my $seq_id = $seq_obj->display_id();
    my $species_id = ( split(m/\|/, $seq_id) )[0];
    $seq_hash{$seq_id} = ();
    $species_hash{$species_id} = ();
  }
}
print STDERR "[INFO] Read in ".commify(scalar(keys %seq_hash))." sequences from ".commify(scalar(@fastas))." files\n";

## dump all seqs in unassigned for later assessment
my %unassigned_seqs = %seq_hash;
my $last_og_id;

## open outfile
open (my $OUT, ">$outprefix.reformatted.txt");

## open Orthogroups
open (my $GROUPS, $orthogroups_file) or die $!;
while (my $line = <$GROUPS>) {
  chomp $line;
  next if $. == 1; ## skip first line
  my @a = split(/\s+/, $line);
  my @b = @a[3..$#a]; ## seqs begin in column 4
  my @c = map{s/,//g; $_} @b;
  print $OUT "$a[1]: " . join (" ", @c) . "\n"; ## NB some OG ids will be replicated in the output file
  print STDERR "\r[INFO] Working on HOG \#$.: $a[0]"; $|=1;

  ## delete from %unassigned_seqs
  ## only genes not found in OGs remain in this hash
  foreach (@c) {
    delete $unassigned_seqs{$_};
  }
  $last_og_id = $a[1];
}
close $GROUPS;

## add unassigned genes to end of OUT file
unless ( $no_unassigned ) {
  print STDERR "[INFO] Adding unassigned genes to outfile\n";
  $last_og_id =~ s/OG//; ## get the number

  foreach (nsort keys %unassigned_seqs) {
    print $OUT "OG" . $last_og_id . ": $_\n";
    $last_og_id++;
  }
}
close $OUT;

print STDERR "[INFO] Printing matrix to file '$outprefix.reformatted.txt'\n";
print STDERR "[INFO] Finished on ".`date`."\n";

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
