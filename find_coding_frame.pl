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
my $aa_fh = Bio::SeqIO -> new ( -file => $aa_file, -format => "fasta" );
while ( my $seq_obj = $aa_fh -> next_seq() ) {
  $prot_hash{$seq_obj->display_id()} = $seq_obj;
}
print STDERR "[INFO] Got ".scalar(keys %prot_hash)." protein seqs from '$aa_file'\n";

## get nuc seqs
my %transcripts_hash;
my $transcripts_fh = Bio::SeqIO -> new ( -file => $dna_file, -format => "fasta" );
while ( my $seq_obj = $transcripts_fh -> next_seq() ) {
  $transcripts_hash{$seq_obj->display_id()} = $seq_obj;
}
print STDERR "[INFO] Got ".scalar(keys %transcripts_hash)." protein seqs from '$dna_file'\n";

## trimmed transcripts
my %stats_hash;
my %results_hash;
my ($unchanged,$fr0,$fr1,$fr2) = (0,0,0,0);
# open (my $LOG, ">find_coding.log") or die "$!\n";

## cycle thru gene ids
foreach my $gid (nsort keys %prot_hash) {
  my $pseq_obj = $prot_hash{$gid};
  my $tseq_obj = $transcripts_hash{$gid};
  ## translate frame 0 and remove terminator '*'
  (my $tseq_translation_fr0 = $tseq_obj->translate( -frame => 0 )->seq()) =~ s/\*$//;

  if ( $pseq_obj->seq() ne $tseq_translation_fr0 ) {
    ## get alternative coding frames
    my $tseq_translation_fr1 = $tseq_obj->translate( -frame => 1 )->seq();
    $tseq_translation_fr1 =~ s/\*$//; ## remove terminator '*'
    my $tseq_translation_fr2 = $tseq_obj->translate( -frame => 2 )->seq();
    $tseq_translation_fr2 =~ s/\*$//; ## remove terminator '*'

    ## check if any match exactly
    my ($m0,$m1,$m2) = ('','','');
    if ( $pseq_obj->seq() eq $tseq_translation_fr1 ) {
      ## correct frame is +1
      ## trim 1 bp from start, and N from end to ensure % 3 == 0
      my $trimmed_seq_string = substr($tseq_obj->seq(), 1, (($tseq_obj->length-1) - (($tseq_obj->length-1) % 3)));
      # my $trimmed_seq_obj = Bio::PrimarySeq -> new ( -id => $gid, -seq => $trimmed_seq_string );
      $results_hash{$gid} = $trimmed_seq_string;

      print $tseq_obj->length() % 3 . "\n";
      print length($trimmed_seq) % 3 . "\n";
      print "#aa in prot: ".$pseq_obj->length()."\n";
      print "#codons in trimmed: ".length($trimmed_seq)/3."\n";
      $m1 = "<==";
      $fr1++;

    } elsif ( $pseq_obj->seq() eq $tseq_translation_fr2 ) {
      ## correct frame is +2
      ## trim 2 bp from start, and N from end to ensure % 3 == 0
      my $trimmed_seq_string = substr($tseq_obj->seq(), 2, (($tseq_obj->length-2) - (($tseq_obj->length-2) % 3)));
      # my $trimmed_seq_obj = Bio::PrimarySeq -> new ( -id => $gid, -seq => $trimmed_seq_string );
      $results_hash{$gid} = $trimmed_seq_string;

      print $tseq_obj->length() % 3 . "\n";
      print length($trimmed_seq) % 3 . "\n";
      print "#aa in prot: ".$pseq_obj->length()."\n";
      print "#codons in trimmed: ".length($trimmed_seq)/3."\n";
      $m2 = "<==";
      $fr2++;

    } else {
      ## leave as frame 0
      $results_hash{$gid} = $tseq_obj->seq();
      $m0 = "<==";
      $fr0++;
    }

    print $gid . "\tPRED\t\t" . $pseq_obj->seq() . "\n";
    print $gid . "\tFRA0\t$m0\t" . $tseq_translation_fr0 . "\n";
    print $gid . "\tFRA1\t$m1\t" . $tseq_translation_fr1 . "\n";
    print $gid . "\tFRA2\t$m2\t" . $tseq_translation_fr2 . "\n";
    print "\n";

  } else {
    ## translation is good
    $results_hash{$gid} = $tseq_obj->seq();
    $unchanged++;
  }
}

# foreach my $gid (nsort )