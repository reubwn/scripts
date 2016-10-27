#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use List::Util qw/sum/;

my $usage = "
  Prints assembly info for fasta file(s)
  USAGE: n50.pl <FASTA1> [FASTA2... ]\n
";

die $usage if $ARGV[0] =~ /-h/;

for my $i (0 .. $#ARGV){

	my (@scaff_lengths);
	my ($tot,$num200,$num500,$num1kb,$num10kb,$ass_N50,$ass_N90,$ass_L50,$ass_L90,$nns,$soft,$gc,$nonATGCN) = (0,0,0,0,0,0,0,0,0,0,0,0,0);

	my $in;
	if ($ARGV[0] eq "-"){
		$in  = Bio::SeqIO->new( -fh => \*STDIN, -format => 'fasta' );
	} else {
		$in  = Bio::SeqIO->new( -file => $ARGV[$i], -format => 'fasta' );
	}
	while (	my $seq = $in->next_seq() ){
		push (@scaff_lengths, $seq->length());
		## count NNNs
		my $c = $seq->seq() =~ tr/N//;
		$nns += $c;
		## count G+C (inc soft masked gc)
		my $d = $seq->seq() =~ tr/GCgc//;
		$gc += $d;
		## count soft masked bases
		my $e = $seq->seq() =~ tr/atgc//;
		$soft += $e;
		## count non-ATCGN bases
		my $f = $seq->seq() =~ tr/ATGCNatgcn//c;
		$nonATGCN += $f;
	}

	my $span = sum @scaff_lengths;

	my @sorted_scaffs = sort {$b <=> $a} @scaff_lengths;

	foreach my $len (@sorted_scaffs){
		$num200++ if $len > 200;
		$num500++ if $len > 500;
		$num1kb++ if $len > 1000;
		$num10kb++ if $len > 10000;
	}

	foreach my $len (@sorted_scaffs){
		$tot += $len;
		$ass_L50++;
		if ($tot >= ($span / 2)){
			$ass_N50 = $len;
			last;
		}
	}

	$tot = 0;
	foreach my $len (@sorted_scaffs){ 
        	$tot += $len;
        	$ass_L90++;
        	if ($tot >= (($span / 10)*9)){ 
                	$ass_N90 = $len;
                	last;
        	}
	}

	print "
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ARGV[$i]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Total span: ".commify($span)."
Total number sequences: ".commify(scalar(@sorted_scaffs))."
  \>200nt: ".commify($num200)."
  \>500nt: ".commify($num500)."
  \>1kb  : ".commify($num1kb)."
  \>10kb : ".commify($num10kb)."
Smallest : ".commify($sorted_scaffs[-1])."
Longest  : ".commify($sorted_scaffs[0])."
	
\%GC  : ".percentage($gc,$span)."
NNNs : ".commify($nns)." (".percentage($nns,$span).")
soft : ".commify($soft)." (".percentage($soft,$span).")
other: ".commify($nonATGCN)." (".percentage($nonATGCN,$span).")

N50: ".commify($ass_N50)."
num: ".commify($ass_L50)."
N90: ".commify($ass_N90)."
num: ".commify($ass_L90)."
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*other = non-ATGCN characters
*N50   = scaffold length at 50\% genome
*num   = number scaffolds to get to N50
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

sub percentage {
    my $numerator = $_[0];
    my $denominator = $_[1];
    my $places = "\%.2f"; ## default is two decimal places
    if (exists $_[2]){$places = "\%.".$_[2]."f";};
    my $float = (($numerator / $denominator)*100);
    my $rounded = sprintf("$places",$float);
    return "$rounded\%";
}

