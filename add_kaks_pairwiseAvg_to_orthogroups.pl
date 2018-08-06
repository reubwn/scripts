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
my $max_seqs = 100;
my $min_seqs = 2;

GetOptions (
  'i|orthogroups=s' => \$orthogroups,
  'p|proteins=s'     => \$prot_path,
  'c|cds=s'         => \$cds_path,
  'm|max:i'       => \$max_seqs,
  'n|min:i'       => \$min_seqs,
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
    $annot_hash{$F[0]}{hU} = $F[3];
    $annot_hash{$F[0]}{AI} = (sprintf "%.1f", $F[6]); ## round to 1dp
    $annot_hash{$F[0]}{category} = $F[9]; ## key= geneid; val=INGROUP or OUTGROUP
    $annot_hash{$F[0]}{CHS} = $F[10];
    $annot_hash{$F[0]}{tax} = $F[11];
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
GROUP: while (my $line = <$GROUPS>) {
  chomp $line;
  my @a = split (m/\s+/, $line);
  my $og_name = shift @a; $og_name =~ s/://;
  print STDERR "\r[INFO] Working on OG \#$.: $og_name"; $|=1;

  ## skip groups that are too big or too small
  next GROUP if scalar(@a) > $max_seqs;
  next GROUP if scalar(@a) < $min_seqs;

  ## fetch proteins and print to temp file
  my (%protein_seqs, %cds_seqs);
  @protein_seqs{@a} = @protein_hash{@a};
  open (my $PRO, ">$outdir/clustal.pro") or die $!;
  foreach (keys %protein_seqs) {
    my $new_id = join (" ", $_, (join (":", $annot_hash{$_}{hU}, $annot_hash{$_}{category}, $annot_hash{$_}{tax})));
    print $PRO ">$new_id\n" . $protein_seqs{$_}->seq() . "\n";
  }
  close $PRO;

  ## fetch corresponding cds seqs as hash of Bio::Seq objects
  @cds_seqs{@a} = @cds_hash{@a};

  ## sanity check that keys in %protein_seqs are same as %cds_seqs
  my %cmp = map { $_ => 1 } keys %protein_seqs;
  for my $key (keys %cds_seqs) {
    last unless exists $cmp{$key};
    delete $cmp{$key};
  }
  if (%cmp) {
    die "[ERROR] Mismatch between protein and CDS seqids: check fasta headers\n\n";
  }

  ## run alignment
  if (system ("clustalo --infile=$outdir/clustal.pro --outfile=$outdir/prot_clustalo/$og_name.prot_aln.faa --force --threads=$threads") != 0) { die "[ERROR] Problem with clustalo!\n"; }
  ## fetch alignment and backtranslate to nucleotides
  my $get_prot_aln = Bio::AlignIO -> new( -file=>"$outdir/prot_clustalo/$og_name.prot_aln.faa", -format=>"fasta" );
  my $prot_aln_obj = $get_prot_aln -> next_aln();
  my $dna_aln_obj = aa_to_dna_aln($prot_aln_obj, \%cds_seqs);

  ## print to file
  my $n_hgt;
  open (my $DNA, ">$outdir/dna_alns/$og_name.dna_aln.fna");
  foreach my $seq_obj ($dna_aln_obj->each_seq) {
    (my $trim = $seq_obj->display_id()) =~ s/\/*//;
    my $new_id = join (" ", $trim, (join (":", $annot_hash{$trim}{hU}, $annot_hash{$trim}{category}, $annot_hash{$trim}{tax})));
    print $DNA ">$new_id\n" . $seq_obj->seq() . "\n";
    $n_hgt++ if $annot_hash{$trim}{category} eq "OUTGROUP"; ## count number of HGTc in OG
  }
  close $DNA;

  ## get Ka (Dn), Ks (Ds) values
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
}
close $GROUPS;
close $OUT;

system ("rm $outdir/clustal*");
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
