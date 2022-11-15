#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
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

foreach my $gid (keys %prot_hash) {
  my $pseq_obj = $prot_hash{$gid};
  my $dseq_obj = $transcripts_hash{$gid};
  my $dseq_translation_obj = $dseq_obj -> translate();
  print $gid . "\t" . $dseq_translation_obj->seq() . "\n";
}
