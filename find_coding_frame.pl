#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Use for trimming fragmented CDSs to the correct reading frame.
  Protein names must match cDNA (transcripts) names exactly.

OPTIONS [*required]
  -a|--aa       *[FILE] : aa sequences (fasta format)
  -d|--dna      *[FILE] : cDNA sequences (fasta format)
  -o|--outprefix [STR]  : outfile prefix ('<INFILE>_inframe.fna')
  -h|--help             : print this message
\n";

my ($aa_file, $dna_file, $help);
my $outprefix = "_inframe.fna";

GetOptions (
  'a|aa=s'      => \$aa_file,
  'd|dna=s'     => \$dna_file,
  'o|outprefix:s'  => \$outprefix,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $aa_file && $dna_file );

## get protein seqs
my %prot_hash;
my $aa_file_fh = Bio::SeqIO -> new ( -file => $aa_file, -format => "fasta" );
while ( my $seq_obj = $aa_file_fh -> next_seq() ) {
  $prot_hash{$seq_obj->display_id()} = $seq_obj;
}
print STDERR "[INFO] Got ".scalar(keys %prot_hash)." protein seqs from '$aa_file'\n";

## get nuc seqs
my %transcripts_hash;
my $dna_file_fh = Bio::SeqIO -> new ( -file => $dna_file, -format => "fasta" );
while ( my $seq_obj = $dna_file_fh -> next_seq() ) {
  $transcripts_hash{$seq_obj->display_id()} = $seq_obj;
}
print STDERR "[INFO] Got ".scalar(keys %transcripts_hash)." protein seqs from '$dna_file'\n";

foreach my $gid (nsort keys %prot_hash) {
  my $pseq_obj = $prot_hash{$gid};
  my $dseq_obj = $transcripts_hash{$gid};
  my $dseq_translation_0 = $dseq_obj->translate( -frame => 0 )->seq();
  $dseq_translation_0 =~ s/\*$//; ## remove terminator '*'
  # print $gid . "\t" . $pseq_obj->seq() . "\n";
  # print $gid . "\t" . $dseq_translation . "\n";
  # print "\n";

  if ( $pseq_obj->seq() ne $dseq_translation_0 ) {
    ## get alternative coding frames
    my $dseq_translation_1 = $dseq_obj->translate( -frame => 1 )->seq();
    $dseq_translation_1 =~ s/\*$//; ## remove terminator '*'
    my $dseq_translation_2 = $dseq_obj->translate( -frame => 2 )->seq();
    $dseq_translation_2 =~ s/\*$//; ## remove terminator '*'

    ## check if any match exactly
    my ($m0,$m1,$m2) = ('','','');
    my $substring;
    my $modulo;
    if ( $pseq_obj->seq() eq $dseq_translation_1 ) {
      print $dseq_obj->length() % 3 . "\n";
      $substring = substr($dseq_obj->seq(), 1, (($dseq_obj->length-1) - ($dseq_obj->length % 3)));
      print length($substring) % 3 . "\n";
      $m1 = "<==";

    } elsif ( $pseq_obj->seq() eq $dseq_translation_2 ) {
      print $dseq_obj->length() % 3 . "\n";
      $substring = substr($dseq_obj->seq(), 2, (($dseq_obj->length-2) - ($dseq_obj->length % 3)));
      print length($substring) % 3 . "\n";
      $m2 = "<==";

    } else {
      $m0 = "<==";
    }

    print $gid . "\tPRED\t\t" . $pseq_obj->seq() . "\n";
    print $gid . "\tFRA0\t$m0\t" . $dseq_translation_0 . "\n";
    print $gid . "\tFRA1\t$m1\t" . $dseq_translation_1 . "\n";
    print $gid . "\tFRA2\t$m2\t" . $dseq_translation_2 . "\n";
    print "\n";

  }
}
