#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;
use Bio::AlignIO;
use File::Basename;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Data::Dumper;

my $usage = "
SYNOPSIS
  Maps dna sequences onto aa alignments, for Ka/Ks, PAML etc.
  Fasta headers must correspond between aa and dna sequences, otherwise the entry in the aa alignment will be skipped.

OPTIONS [*required]
  -a|--aa    *[PATH] : path to dir of aa alignments (fasta format)
  -d|--dna   *[PATH] : path to dir of unaligned dna sequences (fasta format)
  -e|--ext    [STR]  : filename extension to glob from dna path ('fasta')
  -m|--max    [INT]  : maximum number of seqs in aa alignment, skips if > (100)
  -n|--min    [INT]  : minimum number of seqs in aa alignment (2)
  -o|--outdir [DIR]  : base dirname to write stuff ('aa_to_dna_aln_results/')
  -x|--overwrite     : overwrite outdir if it already exists
  -h|--help          : print this message
\n";

my ($aa_path, $dna_path, $outdir, $overwrite, $help);
my $extension = "fasta";
my $max_seqs = 100;
my $min_seqs = 2;

GetOptions (
  'a|aa=s'      => \$aa_path,
  'd|dna=s'     => \$dna_path,
  'e|ext:s'     => \$extension,
  'm|max:i'     => \$max_seqs,
  'n|min:i'     => \$min_seqs,
  'o|outdir:s'  => \$outdir,
  'x|overwrite' => \$overwrite,
  'h|help'      => \$help
);

die $usage if ( $help );
die $usage unless ( $aa_path && $dna_path );

## parse CDSs
my %cds_hash;
my @dna_files = glob ("$dna_path/*.$extension");
print STDERR "[INFO] Reading files from $dna_path/\*.$extension...\n";
foreach my $dna_file (@dna_files) {
  my $in = Bio::SeqIO->new( -file => $dna_file, -format => 'fasta' );
  while (my $seq = $in->next_seq() ) {
    $cds_hash{$seq->display_id()} = $seq->seq();
  }
}
if ( scalar(keys %cds_hash) == 0 ) {
  die "[ERROR] No sequences found in $dna_path!\n";
} else {
  print STDERR "[INFO] Fetched ".commify(scalar(keys %cds_hash))." CDS seqs from ".commify(scalar(@dna_files))." files in $dna_path\n";
}

## cycle thru alignments
my @aln_files = glob ("$aa_path/*.fa");
print STDERR "[INFO] Reading files from $aa_path/\*.fa...\n";
foreach my $aln_file (@aln_files) {
  ## fetch alignment and backtranslate to nucleotides
  my $get_prot_aln = Bio::AlignIO -> new( -file => $aln_file, -format => 'fasta' );
  my $prot_aln = $get_prot_aln -> next_aln();
  my %cds_seqs;
  foreach my $seq ( $prot_aln->each_seq() ) {
    $cds_seqs{$seq->display_id()} = Bio::Seq->new( -display_id => $seq->display_id(), -seq => $cds_hash{$seq->display_id()} );
  }
  my $dna_aln = aa_to_dna_aln($prot_aln, \%cds_seqs);
  my $dna_aln_filename = (basename ($aln_file, ".fa")) . "_dna.fa";
  my $write_dna_aln = Bio::AlignIO -> new( -file => ">$dna_aln_filename", -format => 'fasta' );
  $write_dna_aln -> write_aln($dna_aln);
}

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
