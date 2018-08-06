#!/usr/bin/env perl

## reubwn Aug 2018

use strict;
use warnings;

use Getopt::Long;
use Term::ANSIColor;
use Sort::Naturally;
use File::Path 'rmtree';
use Algorithm::Combinatorics qw(combinations);

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq::EncodedSeq;
use Bio::Align::DNAStatistics;
use Bio::Align::Utilities qw(aa_to_dna_aln);

my $usage = "
SYNOPSIS
  Calculates average pairwise Ka and Ks for all pairs of genes from an OrthoFinder Orthogroups.txt file.
  Uses MAFFT to generate alignments, make sure MAFFT is discoverable in \$PATH.

OPTIONS
  -i|--orthogroups [FILE] : Orthogroups.txt file from OrthoFinder
  -p|--protein     [FILE] : fasta file of protein sequences
  -c|--cds         [FILE] : fasta file of corresponding CDSs (nucleotide)
  -t|--threads     [INT]  : number of clustalo aligning threads
  -a|--annot       [FILE] : annotate sequences with results from HGT analysis
  -o|--outfile     [STR]  : output filename (default: 'inputfilename.kaks.txt')
  -d|--outdir      [DIR]  : base dirname to write stuff (default: 'inputfilename_kaks')
  -h|--help               : print this message
\n";

my ($orthogroups, $prot_path, $cds_path, $annot, $outdir, $outfile, $overwrite, $help);
my $threads = 1;

GetOptions (
  'i|orthogroups=s' => \$orthogroups,
  'p|proteins=s'     => \$prot_path,
  'c|cds=s'         => \$cds_path,
  't|threads:i'   => \$threads,
  'a|annot:s'       => \$annot,
  'o|outfile:s'    => \$outfile,
  'd|outdir:s'     => \$outdir,
  'x|overwrite'   => \$overwrite,
  'h|help'          => \$help
);

die $usage if $help;
die $usage unless ($orthogroups && $prot_path && $cds_path);

## outfiles
$outdir = $orthogroups."_kaks" unless ($outdir);
$outfile = $orthogroups.".kaks.txt" unless ($outfile);

## make $outdir
if (-e $outdir && -d $outdir) {
  if ($overwrite) {
    rmtree([ "$outdir" ]);
    mkdir $outdir;
    mkdir "$outdir/prot_clustalo";
    mkdir "$outdir/dna_alns";
  } else {
    die "[ERROR] Dir $outdir already exists. Specify '-x' to overwrite it.\n";
  }
} else {
  mkdir $outdir;
  mkdir "$outdir/prot_clustalo";
  mkdir "$outdir/dna_alns";
}

## parse proteins and CDSs
my (%protein_hash, %cds_hash);
my @prot_fastas = @{ get_fasta($prot_path) };
print STDERR "[INFO] Reading protein sequences from:\n";
foreach (@prot_fastas){
  print STDERR "-->" . colored($_, 'white on_blue') . "\n";
  my $in = Bio::SeqIO->new ( -file => $_, -format => "fasta" );
  while ( my $seq_obj = $in->next_seq() ){
    $protein_hash{($seq_obj->display_id())} = $seq_obj;
  }
}
print STDERR "[INFO] Read in ".commify(scalar(keys %protein_hash))." protein sequences\n";
my @cds_fastas = @{ get_fasta($cds_path) };
print STDERR "[INFO] Reading CDS sequences from:\n";
foreach (@cds_fastas){
  print STDERR "-->" . colored($_, 'white on_blue') . "\n";
  my $in = Bio::SeqIO->new ( -file => $_, -format => "fasta" );
  while ( my $seq_obj = $in->next_seq() ){
    $cds_hash{($seq_obj->display_id())} = $seq_obj;
  }
}
print STDERR "[INFO] Read in ".commify(scalar(keys %cds_hash))." CDS sequences\n\n";

## check there are some sequences
if ((scalar(keys %protein_hash) == 0) || (scalar(keys %cds_hash) == 0)) {
  die "[ERROR] No sequences found in $prot_path or $cds_path!\n";
}

## parse $annot if present
my %annot_hash;
if ($annot) {
  print STDERR "[INFO] Collecting annotations from " . colored($annot, 'white on_blue') . "\n";
  open (my $ANNOT, $annot) or die $!;
  while (my $line = <$ANNOT>) {
    chomp $line;
    my @F = split (m/\s+/, $line);
    # $annot_hash{$F[0]}{hU} = $F[3];
    # $annot_hash{$F[0]}{AI} = $F[6];
    $annot_hash{$F[0]} = $F[9]; ## key= geneid; val=INGROUP or OUTGROUP
    # $annot_hash{$F[0]}{CHS} = $F[10];
    # $annot_hash{$F[0]}{tax} = $F[11];
  }
  close $ANNOT;
  print STDERR "[INFO] Collected annotations for ".commify(scalar(keys %annot_hash))." genes\n";
}

###############
## MAIN LOOP ##
###############

open (my $OUT, ">$outdir/$outfile") or die $!;
if ($annot) {
  print $OUT join ("\t", "NAME", "NUM_SEQS", "NUM_HGT", "PROP_HGT", "KA", "KS", "KA_VAR", "KS_VAR", "Z_SCORE") . "\n";
} else {
  print $OUT join ("\t", "NAME", "NUM_SEQS", "KA", "KS", "KA_VAR", "KS_VAR", "Z_SCORE") . "\n";
}

## open groups file
open (my $GROUPS, $orthogroups) or die $!;
while (my $line = <$GROUPS>) {
  chomp $line;
  my @a = split (m/\s+/, $line);
  my $og_name = shift @a; $og_name =~ s/://;
  print STDERR "\r[INFO] Working on OG \#$.: $og_name"; $|=1;

  ## fetch proteins and print to temp file
  my (%protein_seqs, %cds_seqs);
  @protein_seqs{@a} = @protein_hash{@a};
  open (my $PRO, ">clustal.pro") or die $!;
  foreach (keys %protein_seqs) {
    print $PRO ">$_\n" . $protein_seqs{$_}->seq() . "\n";
  }
  close $PRO;

  ## get % HGT genes in OG
  my $n_hgt;
  if ($annot) {
    my $string = join (" ", @annot_hash{keys %protein_seqs});
    print STDERR ": $string\n";
    $n_hgt = () = $string =~ m/OUTGROUP/g;
  }

  ## fetch corresponding cds seqs as hash of Bio::Seq objects
  @cds_seqs{@a} = @cds_hash{@a};
  #foreach (keys %cds_seqs) { print "$_\n$cds_seqs{$_}\n" };

  ## sanity check that keys in %protein_seqs are same as
  my %cmp = map { $_ => 1 } keys %protein_seqs;
  for my $key (keys %cds_seqs) {
    last unless exists $cmp{$key};
    delete $cmp{$key};
  }
  if (%cmp) {
    die "[ERROR] Mismatch between protein and CDS seqids: check fasta headers\n\n";
  }

  ## run alignment
  print STDERR "[INFO] Running Clustalo: ".`date`."\n";
  if (system ("clustalo --infile=clustal.pro --outfile=$outdir/prot_clustalo/$og_name.prot_aln.faa --force --threads=$threads") != 0) { die "[ERROR] Problem with clustalo!\n"; }
  ## fetch alignment and backtranslate to nucleotides
  my $get_prot_aln = Bio::AlignIO -> new( -file=>"$outdir/prot_clustalo/$og_name.prot_aln.faa", -format=>"fasta" );
  my $write_dna_aln = Bio::AlignIO -> new( -file=>">$outdir/dna_alns/$og_name.dna_aln.fna", -format=>"fasta" );
  my $prot_aln_obj = $get_prot_aln -> next_aln();
  my $dna_aln_obj = aa_to_dna_aln($prot_aln_obj, \%cds_seqs);
  $write_dna_aln -> write_aln($dna_aln_obj);
  print STDERR "[INFO] Finished Clustalo: ".`date`."\n";

  ## get Ka (Dn), Ks (Ds) values
  print STDERR "[INFO] Getting Ka/Ks: ".`date`."\n";
  eval {
    my $stats = Bio::Align::DNAStatistics->new();
    my $result = $stats->calc_average_KaKs($dna_aln_obj, 1000);
    my ($D_n, $D_s, $D_n_var, $D_s_var, $z_score);
    for (sort keys %{$result}) {
      next if /Seq/;
      if($_ eq "D_n"){$D_n = $result->{$_}};
      if($_ eq "D_s"){$D_s = $result->{$_}};
      if($_ eq "D_n_var"){$D_n_var = $result->{$_};}
      if($_ eq "D_s_var"){$D_s_var = $result->{$_};}
      if($_ eq "z_score"){$z_score = $result->{$_};}
    }
    $D_n = -2 unless ($D_n); ## default values
    $D_s = -2 unless ($D_s);

    ## print to file
    if ($annot) {
      print $OUT join ("\t", $og_name, scalar(keys %cds_seqs), $n_hgt, ($n_hgt/scalar(keys %cds_seqs)), $D_n, $D_s, $D_n_var, $D_s_var, $z_score) . "\n";
    } else {
      print $OUT join ("\t", $og_name, scalar(keys %cds_seqs), $D_n, $D_s, $D_n_var, $D_s_var, $z_score) . "\n";
    }
  };
  print STDERR "[INFO] Finished: ".`date`."\n";
}
close $GROUPS;
close $OUT;

system ("rm clustal*");
print STDERR "\n[INFO] Finished on ".`date`."\n";

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

sub get_fasta {
  my $path = $_[0];
  my @files = glob("$path/*\.f*a");
  if (scalar(@files) == 0) {
    die "[ERROR] Nothing with extension fasta|faa|fna|fa found in $path\n";
  }
  return \@files;
}

__END__
