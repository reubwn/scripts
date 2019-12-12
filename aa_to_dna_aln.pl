#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Data::Dumper;

my $usage = "
SYNOPSIS
  Maps dna sequences onto aa alignments, for Ka/Ks, PAML etc.
  Fasta headers must correspond between aa and dna sequences, otherwise the entry in the aa alignment will be skipped.

OPTIONS [*required]
  -a|--aa    *[PATH] : path to dir of aa alignments (fasta format)
  -d|--dna   *[PATH] : path to dir of unaligned dna sequences (fasta format)
  -m|--max    [INT]  : maximum number of seqs in aa alignment, skips if > (100)
  -n|--min    [INT]  : minimum number of seqs in aa alignment (2)
  -o|--outdir [DIR]  : base dirname to write stuff ('aa_to_dna_aln_results/')
  -x|--overwrite     : overwrite outdir if it already exists
  -h|--help          : print this message
\n";

my ($aa_path, $dna_path, $outdir, $overwrite, $help);
my $max_seqs = 100;
my $min_seqs = 2;

GetOptions (
  'a|aa=s'      => \$aa_path,
  'd|dna=s'     => \$dna_path,
  'm|max:i'     => \$max_seqs,
  'n|min:i'     => \$min_seqs,
  'd|outdir:s'  => \$outdir,
  'x|overwrite' => \$overwrite,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $aa_path && $dna_path );

## parse CDSs
my %cds_hash;
my @files = glob ("$dna_path/*.fna");
foreach my $file (@files) {
  my $in = Bio::SeqIO->new( -file => $file, -format => 'fasta' );
  while (my $seq = $in->next_seq() ) {
    $cds_hash{$seq->display_id()} = $seq->seq();
  }
}
print STDERR "[INFO] Fetched ".commify(scalar(keys %cds_hash))." CDS seqs from ".commify(scalar(@files))." files in $dna_path\n";
die "[ERROR] No sequences found in $dna_path!\n" if ( scalar(keys %cds_hash) == 0 );
