#!/usr/bin/env perl

## Author: reubwn Aug 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "
bam2tab.pl
##########

Converts SAM/BAM file to TAB format used in SSPACE.
Optionally filters reads based on mapping distance from edge of contigs, to remove potential PE contaminant reads from MP datasets.
Calculates mapping edge distances for pairs of reads, excludes pairs with edge distance < threshold (and optional edge distance > threshold):

  contig1                      R2--->
  5++++++====++++++ /break/ +++++====+++++++++3
         <---R1                         contig2

  The regions marked +++ must be > min_edge_distance && < max_edge_distance

NOTE: requires samtools in \$PATH and SAM/BAM sorted by readname (-n option in samtools sort).
USAGE: bam2tab.pl -i <bam_file> [-d mapping_distance] [-s] [-o output_table.txt]

OPTIONS:
  -i|--in                [FILE]  : SAM/BAM file [required]
  -f|--fasta             [FILE]  : fasta file of contigs [required]
  -d|--min_edge_distance [INT]   : discard reads with edge distance < INT [default: 500]
  -m|--max_edge_distance [INT]   : discard reads with edge distance > INT [default: 10,000]
  -s|--same_contig               : discard reads mapping to the same contig [default: no]
  -o|--out               [STR]   : output prefix [default: mapping_table]
  -t|--outtype           [STR]   : Output format: 'table' => table [default], 'sam' => sam format
  -h|--help                      : prints this help message

EXAMPLES:
  (1) Discard reads that map to within 1 kb of contig edges and that map to the same contig:
      >> bam2tab.pl -i mapping.bam -d 1000 -s -o lib1
\n";

## params with defaults
my $min_edge_distance = 500;
my $max_edge_distance = 10000;
my $output = "mapping_table";
my $outtype = "table";

## other args
my ($sam_file,$fasta,$same,$prefix,$help);

GetOptions (
  'in|i=s'                 => \$sam_file,
  'fasta|f=s'              => \$fasta,
  'min_edge_distance|d:i'  => \$min_edge_distance,
  'max_edge_distance|m:i'  => \$max_edge_distance,
  'same_contig|s'          => \$same,
  'out|o:s'                => \$prefix,
  'outtype|t:s'            => \$outtype,
  'help|h'                 => \$help,
);

die $usage if $help;
die $usage unless ($sam_file && $fasta);

## check samtools in $PATH
if ( system("samtools view &>/dev/null" ) == -1){
  die "[ERROR] Samtools error: is samtools in \$PATH?";
}

## set prefix if exists
$output = "$prefix.$output" if $prefix;

## get sequence lengths from fasta
my %lengths;
print "\nFetching contig lengths from $fasta...\n";
my $in = Bio::SeqIO->new ( -file => $fasta, -format => "fasta" );
while ( my $seq_obj = $in->next_seq() ){
  $lengths{($seq_obj->display_id())} = ($seq_obj->length());
}
print "Found ".commify(scalar(keys %lengths))." contigs\n\n";

## open SAM/BAM file, only include proper pairs
my $SAM;
if ($outtype =~ m/sam/i){
  $output = "$output.sam";
  open ($SAM, "samtools view -h -f1 -F3340 $sam_file |") or die $!;
  print "Output set to SAM\n";
} else {
  $output = "$output.tab";
  open ($SAM, "samtools view -f1 -F3340 $sam_file |") or die $!;
  print "Output set to TAB\n";
}
print "MIN mapping distance set to ".commify($min_edge_distance)." nt\n";
print "MAX mapping distance set to ".commify($max_edge_distance)." nt\n\n";

my ($processed,$printed) = (0,0);
open (my $OUT, ">$output") or die $!;

while (<$SAM>){

  ## skip headers unless print to SAM
  if ($_ =~ m/^\@/){
    if ($outtype =~ m/sam/i){
      print $OUT $_;
      next;
    } else {
      next;
    }
  }

  my @read1 = split (/\t/, $_); ## first read
  my @read2 = split (/\t/, <$SAM>); ## second read

  ## skip reads mapping to the same contig
  if ( ($same) && ($read1[6] eq "\=") ){
    $processed++; ## still count them
    next;
  }

  ############################ calculate distances from contig edge

  #my ($fragment_leftmost,$fragment_rightmost,$distance_from_left,$distance_from_right);
#  my ($leftPos1,$rightPos1,$contig1, $leftPos2,$rightPos2,$contig2);
#  my ($leftEgde1,$rightEdge1,$leftEdge2,$rightEdge2) = (0,0,0,0);
  my ($alnLength1,$alnLength2) = (0,0);

  ## get read length from CIGAR
  my ($cigar1,$cigar2) = ($read1[5],$read2[5]);
  $cigar1 =~ s/(\d+)[MX=DN]/$alnLength1+=$1/eg;
  $cigar2 =~ s/(\d+)[MX=DN]/$alnLength2+=$1/eg;

  ## get left/right mapping coords for both reads
  ## NB not concerned with read orientation, ie $leftPos1 may be read1 start or end
  my $leftPos1 = $read1[3];
  my $rightPos1 = $leftPos1 + $alnLength1;
  my $contig1 = $read1[2];

  my $leftPos2 = $read2[3];
  my $rightPos2 = $leftPos2 + $alnLength2;
  my $contig2 = $read2[2];

  my $leftEdge_fromContigStart1 = $leftPos1;
  my $leftEdge_fromContigEnd1 = $lengths{$contig1} - $leftPos1;
  my $rightEdge_fromContigStart1 = $leftPos1 + $alnLength1;
  my $rightEdge_fromContigEnd1 = $lengths{$contig1} - $rightPos1;

  my $leftEdge_fromContigStart2 = $leftPos2;
  my $leftEdge_fromContigEnd2 = $lengths{$contig2} - $leftPos2;
  my $rightEdge_fromContigStart2 = $leftPos2 + $alnLength2;
  my $rightEdge_fromContigEnd2 = $lengths{$contig2} - $rightPos2;


  # ## read1 is on reverse strand, mate is on forward strand
  # if ($read1[1]&16) {
  #   my ($aln_length_1,$aln_length_2) = (0,0);
  #   my ($cigar1,$cigar2) = ($read1[5],$read2[5]);
  #   $cigar1 =~ s/(\d+)[MX=DN]/$aln_length_1+=$1/eg;
  #   $cigar2 =~ s/(\d+)[MX=DN]/$aln_length_2+=$1/eg;
  #
  #   ##                                       R2
  #   ## contig1                               ---->
  #   ## 5'======+++++++++++++ /break/ ++++++++=========3' F strand
  #   ## 3'======+++++++++++++ /break/ ++++++++=========5' R strand
  #   ##    <----                      contig2
  #   ##    R1
  #   ##
  #   ## the regions marked +++ must be > $min_edge_distance && < $max_edge_distance
  #
  #   ## get start/end coords for both reads
  #   $start1 = $read1[3] + $aln_length_1;
  #   $end1 = $read1[3];
  #   $contig1 = $read1[2];
  #   $start2 = $read2[3];
  #   $end2 = $read2[3] + $aln_length_2;
  #   $contig2 = $read2[2];
  #
  #   $edge1 = $lengths{$contig1} - $start1;
  #   $edge2 = $start2;
  #
  # ## read2 is on reverse strand, mate is on forward strand
  # } elsif ($read2[1]&16){
  #
  # ## calculate aln length from CIGAR
  # my ($aln_length_1,$aln_length_2) = (0,0);
  # my ($cigar1,$cigar2) = ($read1[5],$read2[5]);
  # $cigar1 =~ s/(\d+)[MX=DN]/$aln_length_1+=$1/eg;
  # $cigar2 =~ s/(\d+)[MX=DN]/$aln_length_2+=$1/eg;
  #
  # ##                                       R1
  # ## contig1                               ---->
  # ## 5'======+++++++++++++ /break/ ++++++++=========3' F strand
  # ## 3'======+++++++++++++ /break/ ++++++++=========5' R strand
  # ##    <----                      contig2
  # ##    R2
  # ##
  # ## the regions marked +++ must be > $min_edge_distance && < $max_edge_distance
  #
  # ## get start/end coords for both reads
  # $start1 = $read1[3];
  # $end1 = $read1[3] + $aln_length_1;
  # $contig1 = $read1[2];
  # $start2 = $read2[3] + $aln_length_2;
  # $end2 = $read2[3];
  # $contig2 = $read2[2];
  #
  # $edge1 = $start1;
  # $edge2 = $lengths{$contig2} - $start2;
  #
  # }

  ## print to $OUT if all edge distances > $min_edge_distance
  if ( ($leftEdge_fromContigStart1 > $min_edge_distance) &&
       ($leftEdge_fromContigEnd1) > $min_edge_distance &&
       ($rightEdge_fromContigStart1) > $min_edge_distance &&
       ($rightEdge_fromContigEnd1) > $min_edge_distance &&
       ($leftEdge_fromContigStart2) > $min_edge_distance &&
       ($leftEdge_fromContigEnd2) > $min_edge_distance &&
       ($rightEdge_fromContigStart2) > $min_edge_distance &&
       ($rightEdge_fromContigEnd2) > $min_edge_distance ){

    if ($outtype =~ m/sam/i){
      print $OUT join "\t", @read1;
      print $OUT join "\t", @read2;
    } else {
      print $OUT join "\t", $contig1,$leftPos1,$rightPos1,$contig2,$leftPos2,$rightPos2,"\n";
    }
    $printed++;
  }

  ## progress
  $processed++;
  if ($processed % 100000 == 0){
    print "\rProcessed ".commify($processed)." pairs...";
    $| = 1;
  }
}

print "\rProcessed ".commify($processed)." pairs\n";
print "Printed ".commify($printed)." pairs\n\n";

###################### sub-routines

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
