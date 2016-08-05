#!/usr/bin/env perl

## Author: reubwn

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

##
## TODO - set max insert as well? Ie remove pairs with $edge_distance > some_value?
##
##

my $usage = "
bam2tab.pl
Converts SAM/BAM file to TAB format used in SSPACE.
NOTE: requires samtools in \$PATH and SAM/BAM sorted by readname (-n option in samtools sort).

USAGE: bam2tab.pl -i <bam_file> [-d mapping_distance] [-s] [-o output_table.txt]

OPTIONS:
  -i|--in            [FILE]  : SAM/BAM file [required]
  -f|--fasta         [FILE]  : fasta file of contigs [required]
  -d|--edge_distance [INT]   : filter based on minimum mapping distance from contig edge [default: 0]
  -s|--same_contig           : filter reads mapping to the same contig (i.e., are redundant) [default: no]
  -o|--out           [STR]   : output prefix [PREFIX.mapping_table.tab]
  -t|--outtype       [STR]   : Output format: 'table' => table [default], 'sam' => sam format
  -h|--help                  : prints this help message

EXAMPLES:
  (1) Discard reads that map to within 1 kb of contig edges and that map to the same contig:
      >> bam2tab.pl -i mapping.bam -d 1000 -s -p lib1
\n";

## params with defaults
my $edge_distance = 0;
my $output = "mapping_table";
my $outtype = "table";

## other args
my ($sam_file,$fasta,$same,$prefix,$help);

GetOptions (
  'in|i=s'             => \$sam_file,
  'fasta|f=s'          => \$fasta,
  'edge_distance|d:i'  => \$edge_distance,
  'same_contig|s'      => \$same,
  'out|o:s'            => \$prefix,
  'outtype|t:s'        => \$outtype,
  'help|h'             => \$help,
);

die $usage if $help;
die $usage unless ($sam_file && $fasta);

## check samtools in $PATH
if (system("samtools view &>/dev/null")==-1){
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
print "Minimum mapping distance set to ".commify($edge_distance)." nt\n\n";

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

  my @F = split (/\t/, $_); ## first read
  my @R = split (/\t/, <$SAM>); ## second read

  ## skip reads mapping to the same contig
  if ( ($same) && ($F[6] eq "\=") ){
    next;
  }

  ############################ calculate distances from contig edge

  #my ($fragment_leftmost,$fragment_rightmost,$distance_from_left,$distance_from_right);
  my ($start1,$end1,$contig1, $start2,$end2,$contig2);
  my ($edge1,$edge2) = (0,0);

  ## read1 is on forward strand, mate is on reverse strand
  if ($F[1]&32){

    ## calculate aln length from CIGAR
    my ($aln_length_1,$aln_length_2) = (0,0);
    my ($cigar1,$cigar2) = ($F[5],$R[5]);
    $cigar1 =~ s/(\d+)[MX=DN]/$aln_length_1+=$1/eg;
    $cigar2 =~ s/(\d+)[MX=DN]/$aln_length_2+=$1/eg;

    ## contig2
    ## ==========+++++++++++++
    ##    3'<----5'
    ##      read2(R)            contig1
    ##                          ++++++++==================
    ##                                5'---->3'
    ##                                  read1(F)
    ##
    ## the regions marked +++ must be > $edge_distance

    ## get start/end coords for both reads
    $start1 = $F[3];
    $end1 = $F[3] + $aln_length_1;
    $contig1 = $F[2];
    $start2 = $R[3] + $aln_length_2;
    $end2 = $R[3];
    $contig2 = $R[2];

    $edge1 = $start1;
    $edge2 = $lengths{$contig2} - $start2;

  ## read1 is on reverse strand, mate is on forward strand
  } elsif ($F[1]&16) {
    my ($aln_length_1,$aln_length_2) = (0,0);
    my ($cigar1,$cigar2) = ($F[5],$R[5]);
    $cigar1 =~ s/(\d+)[MX=DN]/$aln_length_1+=$1/eg;
    $cigar2 =~ s/(\d+)[MX=DN]/$aln_length_2+=$1/eg;

    ## contig1
    ## ==========+++++++++++++
    ##    3'<----5'
    ##      read1(F)            contig2
    ##                          ++++++++==================
    ##                                5'---->3'
    ##                                     read2(R)
    ##
    ## the regions marked +++ must be > $edge_distance

    ## get start/end coords for both reads
    $start1 = $F[3] + $aln_length_1;
    $end1 = $F[3];
    $contig1 = $F[2];
    $start2 = $R[3];
    $end2 = $R[3] + $aln_length_2;
    $contig2 = $R[2];

    $edge1 = $lengths{$contig1} - $start1;
    $edge2 = $start2;
  }

  ## print to $OUT if both $edge1 && $edge2 > $edge_distance
  if (($edge1 > $edge_distance) && ($edge2 > $edge_distance)){
    if ($outtype =~ m/sam/i){
      print $OUT join "\t", @F;
      print $OUT join "\t", @R;
    } else {
      print $OUT join "\t", $contig1,$start1,$end1,$contig2,$start2,$end2,"\n";
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
