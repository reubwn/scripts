#!/usr/bin/env perl

## reubwn Aug 2018

use strict;
use warnings;

use Getopt::Long;
use Term::ANSIColor;
use Sort::Naturally;
use File::Path qw(rmtree);

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);

my $usage = "
SYNOPSIS
  Aligns orthologous groups from proteins and backtranslate to DNA.
  Uses Clustalo to generate alignments, make sure Clustalo is discoverable in \$PATH.

OPTIONS
  -i|--orthogroups [FILE] : Orthogroups.txt file from OrthoFinder
  -p|--protein     [FILE] : path to protein fasta sequences
  -c|--cds         [FILE] : path to corresponding CDS sequences
  -m|--max         [INT]  : maximum number of seqs in OG, skips if > (100)
  -n|--min         [INT]  : minimum number of seqs in OG (2)
  -t|--threads     [INT]  : number of clustalo aligning threads
  -a|--annot       [FILE] : annotate sequences with results from HGT analysis
  -u|--hU          [INT]  : hU threshold for defining HGTc (>30)
  -s|--CHS         [INT]  : CHS threshold for defining HGTc (>90\%)
  -d|--outdir      [DIR]  : base dirname to write stuff ('prepare/')
  -o|--outfile     [STR]  : base filename for some simple stats ('prepare_stats')
  -x|--overwrite          : overwrite outdir and outfile
  -h|--help               : print this message
\n";

my ($orthogroups, $prot_path, $cds_path, $annot, $outdir, $outfile, $overwrite, $help);
my $threads = 1;
my $max_seqs = 100;
my $min_seqs = 2;
my $hgt_threshold = 30;
my $chs_threshold = 90;

GetOptions (
  'i|orthogroups=s' => \$orthogroups,
  'p|proteins=s'     => \$prot_path,
  'c|cds=s'         => \$cds_path,
  'm|max:i'       => \$max_seqs,
  'n|min:i'       => \$min_seqs,
  't|threads:i'   => \$threads,
  'a|annot:s'       => \$annot,
  'u|hU_threshold:i'  => \$hgt_threshold,
  's|CHS_threshold:i' => \$chs_threshold,
  'o|outfile:s'    => \$outfile,
  'd|outdir:s'     => \$outdir,
  'x|overwrite'   => \$overwrite,
  'h|help'        => \$help
);

die $usage if $help;
die $usage unless ($orthogroups && $prot_path && $cds_path);

## outfiles
$outdir = "prepare" unless ($outdir);
$outfile = "prepare_stats" unless ($outfile);

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
  print STDERR "> " . colored($_, 'white on_blue') . "\n";
  my $in = Bio::SeqIO->new ( -file => $_, -format => "fasta" );
  while ( my $seq_obj = $in->next_seq() ){
    $protein_hash{($seq_obj->display_id())} = $seq_obj;
  }
}
print STDERR "[INFO] Read in ".commify(scalar(keys %protein_hash))." protein sequences\n";
my @cds_fastas = @{ get_fasta($cds_path) };
print STDERR "[INFO] Reading CDS sequences from:\n";
foreach (@cds_fastas){
  print STDERR "> " . colored($_, 'white on_blue') . "\n";
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
  LINE: while (my $line = <$ANNOT>) {
    next LINE if $line =~ m/^\#/;
    chomp $line;
    my @F = split (m/\s+/, $line);
    $annot_hash{$F[0]}{hU} = sprintf("%.1f", $F[3]);
    $annot_hash{$F[0]}{AI} = sprintf("%.1f", $F[6]); ## round to 1dp
    $annot_hash{$F[0]}{category} = $F[9]; ## key= geneid; val=INGROUP or OUTGROUP
    $annot_hash{$F[0]}{CHS} = $F[10];
    $annot_hash{$F[0]}{tax} = $F[11];
  }
  close $ANNOT;
  print STDERR "[INFO] Collected annotations for ".commify(scalar(keys %annot_hash))." genes\n";
}

# G		Glycine		    Gly		P		Proline		    Pro
# A		Alanine		    Ala		V		Valine		    Val
# L		Leucine		    Leu		I		Isoleucine	  Ile
# M		Methionine		Met		C		Cysteine		  Cys
# F		Phenylalanine	Phe		Y		Tyrosine		  Tyr
# W		Tryptophan		Trp		H		Histidine		  His
# K		Lysine        Lys		R		Arginine		  Arg
# Q		Glutamine		  Gln		N		Asparagine	  Asn
# E		Glutamic Acid	Glu		D		Aspartic Acid	Asp
# S		Serine		    Ser		T		Threonine		  Thr
my @acids = qw(G A L M F W K Q E S P V I C Y H R N D T);

open (my $OUTG, ">$outdir/$outfile.groups.txt") or die $!;
open (my $OUTD, ">$outdir/$outfile.dna.txt") or die $!;
open (my $OUTP, ">$outdir/$outfile.protein.txt") or die $!;
if ($annot) {
  print $OUTG join ("\t", "#NAME", "NUM_SEQS", "NUM_HGT", "PROP_HGT") . "\n";
  print $OUTD join ("\t", "#NAME", "OG", "HGTc", "hU", "AI", "BBSUMCAT", "CHS", "TAX", "GC") . "\n";
  print $OUTP join ("\t", "#NAME", "OG", "HGTc", "hU", "AI", "BBSUMCAT", "CHS", "TAX", nsort(@acids)) . "\n";
} else {
  print $OUTG join ("\t", "#NAME", "NUM_SEQS") . "\n";
  print $OUTD join ("\t", "#NAME", "OG", "GC") . "\n";
  print $OUTP join ("\t", "#NAME", "OG", nsort(@acids)) . "\n";
}

###############
## MAIN LOOP ##
###############

## HGTc hash for storing HGTc gene names
my %HGTc_hash;

## open groups file
open (my $GROUPS, $orthogroups) or die $!;

## groups loop
GROUP: while (my $line = <$GROUPS>) {
  chomp $line;
  my @a = split (m/\s+/, $line);
  my $og_name = shift @a; $og_name =~ s/://;
  print STDERR "\r[INFO] Working on OG \#$.: $og_name"; $|=1;

  ## fetch proteins, analyse and print to temp file
  my (%protein_seqs, %cds_seqs);
  @protein_seqs{@a} = @protein_hash{@a};
  open (my $PRO, ">$outdir/clustal.pro") or die $!;
  foreach my $pid (keys %protein_seqs) {
    my $header;
    ## if gene has some annotations in %annot_hash...
    if ($annot_hash{$pid}{category}) {
      ## if gene is HGTc (must have hU > 30 && category eq OUTGROUP with support > 90%)
      if ( ($annot_hash{$pid}{hU} > $hgt_threshold) and ($annot_hash{$pid}{category} eq "OUTGROUP") and ($annot_hash{$pid}{CHS} > $chs_threshold) ) {
        ## gene is HGTc... gets '1'
        ## print details to $OUTP
        print $OUTP join ("\t", $pid, $og_name, "1", $annot_hash{$pid}{hU}, $annot_hash{$pid}{AI}, $annot_hash{$pid}{category}, $annot_hash{$pid}{CHS}, $annot_hash{$pid}{tax});
        ## count residues
        foreach my $res (nsort @acids) {
          my $string = $protein_seqs{$pid}->seq();
          my $count = eval "\$string =~ tr/$res//";
          print $OUTP "\t" . sprintf("%.4f", $count/length($string));
        }
        print $OUTP "\n";
        ## annotate fasta headers
        $header = join (" ", $pid, (join (":", $annot_hash{$pid}{hU}, $annot_hash{$pid}{category}, $annot_hash{$pid}{tax})));
        ## put $pid into %HGTc_hash
        $HGTc_hash{$pid}++;

      } else {
        ## gene is not HGTc... gets '0' but still has some annotations
        ## print details to $OUTP
        print $OUTP join ("\t", $pid, $og_name, "0", $annot_hash{$pid}{hU}, $annot_hash{$pid}{AI}, $annot_hash{$pid}{category}, $annot_hash{$pid}{CHS}, $annot_hash{$pid}{tax});
        ## count residues
        foreach my $res (nsort @acids) {
          my $string = $protein_seqs{$pid}->seq();
          my $count = eval "\$string =~ tr/$res//";
          print $OUTP "\t" . sprintf("%.4f", $count/length($string));
        }
        print $OUTP "\n";
        ## annotate fasta headers
        $header = $pid;
      }
    } else { ## if no annotations present
      ## gene is not HGTc... gets '0' with no annotations
      ## print details to $OUTP
      print $OUTP join ("\t", $pid, $og_name, "0", "NA","NA","NA","NA","NA");
      ## count residues
      foreach my $res (nsort @acids) {
        my $string = $protein_seqs{$pid}->seq();
        my $count = eval "\$string =~ tr/$res//";
        print $OUTP "\t" . sprintf("%.4f", ($count/length($string)));
      }
      print $OUTP "\n";
      ## simple header
      $header = $pid;
    }
    ## print to fasta format
    print $PRO ">$header\n" . $protein_seqs{$pid}->seq() . "\n";
  }
  close $PRO;

  ## fetch corresponding cds seqs as hash of Bio::Seq objects
  @cds_seqs{@a} = @cds_hash{@a};
  foreach my $gid (keys %cds_seqs) {
    ## if gene has annotation...
    if ($annot_hash{$gid}{category}) {
      ## if gene is HGTc (must have hU > 30 && category eq OUTGROUP with support > 90%)
      if ( ($annot_hash{$gid}{hU} > $hgt_threshold) and ($annot_hash{$gid}{category} eq "OUTGROUP") and ($annot_hash{$gid}{CHS} > $chs_threshold) ) {
        ## gene is HGTc... gets '1'
        ## print details to $OUTD
        print $OUTD join ("\t", $gid, $og_name, "1", $annot_hash{$gid}{hU}, $annot_hash{$gid}{AI}, $annot_hash{$gid}{category}, $annot_hash{$gid}{CHS}, $annot_hash{$gid}{tax});
        ## count GC
        my $GC = $cds_seqs{$gid}->seq() =~ tr/GCgc//;
        print $OUTD "\t" . sprintf("%.4f", ($GC/$cds_seqs{$gid}->length())) . "\n";
      } else {
        ## gene is not HGTc... gets '0'
        ## print details to $OUTD
        print $OUTD join ("\t", $gid, $og_name, "0", $annot_hash{$gid}{hU}, $annot_hash{$gid}{AI}, $annot_hash{$gid}{category}, $annot_hash{$gid}{CHS}, $annot_hash{$gid}{tax});
        ## count GC
        my $GC = $cds_seqs{$gid}->seq() =~ tr/GCgc//;
        print $OUTD "\t" . sprintf("%.4f", ($GC/$cds_seqs{$gid}->length())) . "\n";
      }
    } else { ## if no annotations present
    ## gene is not HGTc... gets '0'
    ## print details to $OUTD
    print $OUTD join ("\t", $gid, $og_name, "0", "NA","NA","NA","NA","NA");
    ## count GC
    my $GC = $cds_seqs{$gid}->seq() =~ tr/GCgc//;
    print $OUTD "\t" . sprintf("%.4f", ($GC/$cds_seqs{$gid}->length())) . "\n";
    }
  }

  ## sanity check that keys in %protein_seqs are same as %cds_seqs
  my %cmp = map { $_ => 1 } keys %protein_seqs;
  for my $key (keys %cds_seqs) {
    last unless exists $cmp{$key};
    delete $cmp{$key};
  }
  if (%cmp) {
    die "[ERROR] Mismatch between protein and CDS seqids: check fasta headers\n\n";
  }

  ## print some group stats to file
  if ($annot) {
    print $OUTG join ("\t", $og_name, scalar(@a), scalar(keys %HGTc_hash), (scalar(keys %HGTc_hash)/scalar(@a))) . "\n";
  } else {
    print $OUTG join ("\t", $og_name, scalar(@a)) . "\n";
  }

  ## at this point skip groups that are too big or too small
  next GROUP if scalar(@a) > $max_seqs;
  next GROUP if scalar(@a) < $min_seqs;

  ## run alignment
  if (system ("clustalo --infile=$outdir/clustal.pro --outfile=$outdir/prot_clustalo/$og_name.prot_aln.faa --force --threads=$threads") != 0) { die "[ERROR] Problem with clustalo!\n"; }
  ## fetch alignment and backtranslate to nucleotides
  my $get_prot_aln = Bio::AlignIO -> new( -file => "$outdir/prot_clustalo/$og_name.prot_aln.faa", -format => "fasta" );
  my $prot_aln_obj = $get_prot_aln -> next_aln();
  my $dna_aln_obj = aa_to_dna_aln($prot_aln_obj, \%cds_seqs);

  ## print to file
  open (my $DNA, ">$outdir/dna_alns/$og_name.dna_aln.fna");
  foreach my $seq_obj ($dna_aln_obj->each_seq) {
    (my $temp = $seq_obj->display_id()) =~ s/\/*//; ##trim annoying length suffix
    my $header;
    if ($annot_hash{$temp}{category}) {
      $header = join (" ", $temp, (join (":", $annot_hash{$temp}{hU}, $annot_hash{$temp}{category}, $annot_hash{$temp}{tax})));
    } else {
      $header = $temp;
    }
    print $DNA ">$header\n" . $seq_obj->seq() . "\n";
  }
  close $DNA;
}

close $GROUPS;
close $OUTG;
close $OUTP;
close $OUTD;

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
