#!/usr/bin/env perl

## author: reubwn October 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Term::ANSIColor;
use Sort::Naturally;
use File::Path 'rmtree';

my $usage = "
SYNOPSIS
  Converts Orthogroups.txt to fasta.

OPTIONS
  -i|--groups [FILE] : Orthogroups.txt file [required]
  -p|--path   [FILE] : path to directory of fasta files used to construct Orthogroups.txt [required]
  -d|--outdir [STR]  : output dir name [default: {infile}_seqs]
  -a|--annot  [FILE] : annotate sequences with results from HGT analysis
  -h|--help          : prints this help message
\n";

my ($orthogroups_file,$path,$outdir,$annot,$help);
#my $outdir = "orthofinder_groups2fasta";

GetOptions (
  'i|in=s'     => \$orthogroups_file,
  'p|path=s'   => \$path,
  'd|outdir:s' => \$outdir,
  'a|annot:s'  => \$annot,
  'h|help'     => \$help,
);

die $usage if $help;
die $usage unless ($orthogroups_file && $path);

## get sequences
my %seq_hash;
my @fastas = glob("$path*fasta");
if (scalar(@fastas) == 0) {
  print STDERR "[INFO] Nothing found in $path with *.fasta... will try *.faa\n";
  @fastas = glob("$path*faa");
}
if (scalar(@fastas) == 0) {
  print STDERR "[INFO] Nothing found in $path with *.faa... will try *.fna\n";
  @fastas = glob("$path*fna");
}
if (scalar(@fastas) == 0) {
  print STDERR "[INFO] Nothing found in $path with *.fna... will try *.aa (augustus style)\n";
  @fastas = glob("$$path*aa");
}
if (scalar(@fastas) == 0) {
  die "[ERROR] Still nothing found in $path\nPlease make sure there are protein fastas in $path with extension fasta|faa|fna|aa\n";
}
print STDERR "[INFO] Reading sequences from:\n";
foreach (@fastas){
  print STDERR colored($_, 'white on_blue') . "\n";
  my $in = Bio::SeqIO->new ( -file => $_, -format => "fasta" );
  while ( my $seq_obj = $in->next_seq() ){
    $seq_hash{($seq_obj->display_id())} = ($seq_obj->seq());
  }
}
print STDERR "[INFO] Read in ".commify(scalar(keys %seq_hash))." sequences from ".commify(scalar(@fastas))." files\n\n";

## parse $annot if present
my %annot_hash;
if ($annot) {
  print STDERR "[INFO] Collecting annotations from " . colored($annot, 'white on_blue') . "\n";
  open (my $ANNOT, $annot) or die $!;
  while (my $line = <$ANNOT>) {
    chomp $line;
    my @F = split (m/\s+/, $line);
    $annot_hash{$F[0]}{hU} = $F[3];
    $annot_hash{$F[0]}{AI} = $F[6];
    $annot_hash{$F[0]}{category} = $F[9];
    $annot_hash{$F[0]}{CHS} = $F[10];
    $annot_hash{$F[0]}{tax} = $F[11];
  }
}
close $ANNOT;
print STDERR "[INFO] Collected annotations for ".commify(scalar(keys %annot_hash))." genes\n";

## make $outdir
$outdir = $orthogroups_file."_seqs";
if (-e $outdir && -d $outdir) {
  rmtree([ "$outdir" ]);
  mkdir $outdir;
} else {
  mkdir $outdir;
}

## open Orthogroups
my %orthogroups;
open (my $GROUPS, $orthogroups_file) or die $!;
while (my $line = <$GROUPS>) {
  chomp $line;
  my @a = split(/\:\s+/, $line);
  my @b = split(/\s+/, $a[1]);
  ## open file:
  open (my $OUT, ">$outdir/$a[0].fasta") or die $!;
  foreach (@b) {
    if ($annot) {
      if ($annot_hash{$_}{category} =~ "OUTGROUP") { ## only write annotations for HGTc genes
        print $OUT ">$_ ";
        print $OUT join ("|", $annot_hash{$_}{hU}, $annot_hash{$_}{AI}, $annot_hash{$_}{category}, $annot_hash{$_}{CHS}, $annot_hash{$_}{tax}, "\n")
      } else {
        print $OUT ">$_\n";
      }
      print $OUT "$seq_hash{$_}\n";
    } else {
      print $OUT ">$_\n$seq_hash{$_}\n";
    }
  }
  close $OUT;
  print STDERR "\r[INFO] Working on OG \#$.: $a[0]"; $|=1;
}
close $GROUPS;

print STDERR "\n[INFO] Finished on ".`date`."\n";

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
