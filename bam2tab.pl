#!/usr/bin/env perl

## Author: reubwn

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "
bam2tab.pl
Converts SAM/BAM file to TAB format used in SSPACE.
NOTE: requires samtools in \$PATH and SAM/BAM sorted by readname (-n option in samtools sort).

USAGE: bam2tab.pl -b <bam_file> [-s <sam_file>] [-d mapping_distance] [-o output_table.txt]

OPTIONS:
  -i|--in            [FILE]  : SAM/BAM file [required]
  -f|--fasta         [FILE]  : fasta file of contigs [required]
  -d|--edge_distance [INT]   : filter based on minimum mapping distance from contig edge [default: 0]
  -s|--same_contig           : filter reads mapping to the same contig (i.e., are redundant) [default: no]
  -o|--output        [FILE]  : output file name [default: mapping_table.txt]
  -h|--help                  : prints this help message

EXAMPLES:
  (1) Get reads that map to within 1 kb of contig edges, discard any that map to the same contig:
      >> bam2tab.pl -b mapping.bam -d 1000 -s -o mapping_table.1kb.txt
\n";

## params with defaults
my $edge_distance = 0;
my $output = "mapping_table.txt";

## other args
my ($sam_file,$fasta,$same,$help);

GetOptions (
  'in|i=s'             => \$sam_file,
  'fasta|f=s'          => \$fasta,
  'edge_distance|d:i'  => \$edge_distance,
  'same_contig|s'      => \$same,
  'output|o:s'         => \$output,
  'help|h'             => \$help,
);

die $usage if $help;
die $usage unless ($sam_file && $fasta);

## check samtools in $PATH
if (system("samtools view &>/dev/null")==-1){
  die "[ERROR] Samtools error: is samtools in \$PATH?";
}

## get sequence lengths from fasta
my %lengths;
print "\nFetching contig lengths from $fasta... ";
my $in = Bio::SeqIO->new ( -file => $fasta, -format => "fasta" );
while ( my $seq_obj = $in->next_seq() ){
  $lengths{($seq_obj->display_id())} = ($seq_obj->length());
}
print scalar(keys %lengths)." sequences\n";

open (my $OUT, ">$output") or die $!;

## open SAM/BAM file, skip SAM headers, only include proper pairs
open (my $SAM, "samtools view -f1 -F3340 $sam_file |") or die $!;
while (<$SAM>){
  chomp;
  my @F = split (/\t/, $_); ## first read
  my @R = split (/\t/, <$SAM>); ## second read

  ## skip reads mapping to the same contig
  if ( ($same) && ($F[6] eq "\=") ){
    next;
  }

  ############################ calculate distance from contig edge

  my ($fragment_leftmost,$fragment_rightmost,$distance_from_left,$distance_from_right);

  ## read is on forward strand, mate is on reverse strand
  if ($F[1]&32){

    ## calculate aln length from CIGAR
    my ($aln_length_F,$aln_length_R) = (0,0);
    $F[5] =~ s/(\d+)[MX=DN]/$aln_length_F+=$1/eg;
    $R[5] =~ s/(\d+)[MX=DN]/$aln_length_R+=$1/eg;

    ## calculate outermost fragment positions
    $fragment_leftmost = $F[3];
    $fragment_rightmost = $R[3] + $aln_length_R - 1;

    ## calculate distance from contig edges
    $distance_from_left = $fragment_leftmost;
    $distance_from_right = $lengths{$R[2]} - $fragment_rightmost;

    ## print to $OUT if $distance_from_left AND $distance_from_right > $edge_distance
    if ( ($distance_from_left > $edge_distance) && ($distance_from_right > $edge_distance) ){
      print "$F[2]\t$fragment_leftmost\t".($fragment_leftmost + $aln_length_F - 1)."\t$R[2]\t".($fragment_rightmost - $aln_length_R + 1)."\t$fragment_rightmost\n";
    }

  ## read is on reverse strand, mate is on forward strand
  } elsif ($F[1]&16) {
    my ($aln_length_F,$aln_length_R) = (0,0);
    $F[5] =~ s/(\d+)[MX=DN]/$aln_length_F+=$1/eg;
    $R[5] =~ s/(\d+)[MX=DN]/$aln_length_R+=$1/eg;

    $fragment_leftmost = $R[3];
    $fragment_rightmost = $F[3] + $aln_length_F - 1;

    $distance_from_left = $fragment_leftmost;
    $distance_from_right = $lengths{$F[2]} - $fragment_rightmost;

    if ( ($distance_from_left > $edge_distance) && ($distance_from_right > $edge_distance) ){
      print "$R[2]\t$fragment_leftmost\t".($fragment_leftmost + $aln_length_R - 1)."\t$F[2]\t".($fragment_rightmost - $aln_length_F + 1)."\t$fragment_rightmost\n";
    }
  }
}

__END__
