#!/usr/bin/env perl
use strict;
use warnings;
use Bio::AlignIO;

my $file = shift or die "Usage: $0 aligned_sequences.fasta\n";

my $in = Bio::AlignIO->new(
    -file   => $file,
    -format => 'fasta'
);

my $aln = $in->next_aln;
my @seqs = $aln->each_seq;
my $nseq = scalar @seqs;
my $len  = $aln->length;

die "Need at least two sequences\n" if $nseq < 2;

my $total_differences = 0;
my $total_sites_compared = 0;

for my $i (0 .. $#seqs - 1) {
    for my $j ($i + 1 .. $#seqs) {

        my $s1 = uc $seqs[$i]->seq;
        my $s2 = uc $seqs[$j]->seq;

        my $pair_diffs = 0;
        my $pair_sites = 0;

        for my $pos (0 .. $len - 1) {
            my $a = substr($s1, $pos, 1);
            my $b = substr($s2, $pos, 1);

            # Skip gaps and ambiguous bases
            next unless $a =~ /^[ACGT]$/;
            next unless $b =~ /^[ACGT]$/;

            $pair_sites++;
            $pair_diffs++ if $a ne $b;
        }

        $total_differences += $pair_diffs;
        $total_sites_compared += $pair_sites;
    }
}

die "No comparable sites found\n" if $total_sites_compared == 0;

my $pi = $total_differences / $total_sites_compared;

print "Sequences: $nseq\n";
print "Alignment length: $len\n";
print "Pairwise comparisons: ", $nseq * ($nseq - 1) / 2, "\n";
print "Nucleotide diversity (pi): $pi\n";
