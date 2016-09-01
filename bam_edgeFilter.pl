#!/usr/bin/env perl

## Author: reubwn Aug 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "
bam_edgeFilter.pl
=================

Converts SAM/BAM file to TAB format used in SSPACE.
Optionally filters reads based on mapping distance from edge of contigs, to remove potential PE contaminant reads from MP datasets.
Calculates mapping edge distances for pairs of reads, excludes pairs with edge distance < threshold:

  contig1                      R2--->
  5++++++====++++++ /break/ +++++====+++++++++3
         <---R1                         contig2

  The regions marked +++ must be > min_edge_distance

NOTE: requires samtools in \$PATH and SAM/BAM sorted by readname (-n option in samtools sort).
USAGE: bam_edgeFilter.pl -i <bam_file> -f <reference> [-d INT] [-s] [-o PREFIX]

OPTIONS:
  -i|--in                [FILE] : SAM/BAM file [required]
  -f|--fasta             [FILE] : fasta file of contigs [required]
  -d|--min_edge_distance [INT]  : discard reads with edge distance < INT [default: 500]
  -s|--same_contig              : discard reads mapping to the same contig [default: no]
  -o|--out               [STR]  : output prefix [default: edgeFilter_output]
  -t|--outtype           [STR]  : output format: 'sam' => sam format [default], 'table' => table
  -h|--help                     : prints this help message

EXAMPLES:
  (1) Discard reads that map to within 1 kb of contig edges and that map to the same contig:
      >> bam_edgeFilter.pl -i mapping.bam -f reference.fasta -d 1000 -s -t sam -o lib1.filtered
\n";

## params with defaults
my $min_edge_distance = 500;
my $prefix = "edgeFilter_output";
my $outtype = "sam";

## other args
my ($sam_file,$fasta,$same,$help,$outfile);

GetOptions (
  'in|i=s'                 => \$sam_file,
  'fasta|f=s'              => \$fasta,
  'min_edge_distance|d:i'  => \$min_edge_distance,
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
  $outfile = "$prefix.sam";
  open ($SAM, "samtools view -h -f1 -F3340 $sam_file |") or die $!;
  print "Output set to SAM\n";
} elsif ($outtype =~ m/table/i) {
  $outfile = "$prefix.tab";
  open ($SAM, "samtools view -f1 -F3340 $sam_file |") or die $!;
  print "Output set to TAB\n";
} else {
  die "Outtype not recognised!\n\n";
}
print "Discard reads mapping to same contig set to ";
if ($same) {print "TRUE\n"} else {print "FALSE\n"};
print "MIN mapping distance set to ".commify($min_edge_distance)." nt\n\n";

my ($processed,$printed) = (0,0);
open (my $OUT, ">$outfile") or die $!;

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

  my @read1 = split (/\t/, $_); ## first read in pair
  my @read2 = split (/\t/, <$SAM>); ## second read in pair

  ## skip reads mapping to the same contig
  if ( ($same) && ($read1[6] eq "\=") ){
    $processed++; ## still count them
    next;
  }

  ############################ calculate distances from contig edge

  ## get read length from CIGAR
  my ($alnLength1,$alnLength2) = (0,0);
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

  ## LOGIC:
  ## there are 4 edge distances for each mapped read:
  ##       contig         read
  ## start|++++++++++++<========++++++++++++++++++++++|end
  ##       <----(1)---><------------(2)-------------->
  ##       <--------(3)--------><--------(4)--------->
  ##
  ## thus 8 when both reads are considered
  ## all 8 must be > $min_edge_distance to pass filter

  my $leftEdge_fromContigStart1 = $leftPos1; ## distance (1) above, read1
  my $leftEdge_fromContigEnd1 = $lengths{$contig1} - $leftPos1; ## (2)
  my $rightEdge_fromContigStart1 = $leftPos1 + $alnLength1; ## (3)
  my $rightEdge_fromContigEnd1 = $lengths{$contig1} - $rightPos1; ## (4)

  my $leftEdge_fromContigStart2 = $leftPos2; ## distance (1) above, read2
  my $leftEdge_fromContigEnd2 = $lengths{$contig2} - $leftPos2; ## (2)
  my $rightEdge_fromContigStart2 = $leftPos2 + $alnLength2; ## (3)
  my $rightEdge_fromContigEnd2 = $lengths{$contig2} - $rightPos2; ## (4)

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
  if ($processed % 200000 == 0){
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
