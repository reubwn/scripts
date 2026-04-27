#!/usr/bin/env python3

from Bio import AlignIO
import sys

if len(sys.argv) != 2:
    sys.exit("Usage: python calculate_pi.py alignment.fasta")

alignment = AlignIO.read(sys.argv[1], "fasta")

nseq = len(alignment)
aln_len = alignment.get_alignment_length()

if nseq < 2:
    sys.exit("Need at least two sequences")

total_differences = 0
total_sites_compared = 0

for i in range(nseq - 1):
    for j in range(i + 1, nseq):
        seq1 = str(alignment[i].seq).upper()
        seq2 = str(alignment[j].seq).upper()

        for a, b in zip(seq1, seq2):
            # Skip gaps and ambiguous bases
            if a not in "ACGT" or b not in "ACGT":
                continue

            total_sites_compared += 1

            if a != b:
                total_differences += 1

if total_sites_compared == 0:
    sys.exit("No comparable sites found")

pi = total_differences / total_sites_compared

print(f"Sequences: {nseq}")
print(f"Alignment length: {aln_len}")
print(f"Pairwise comparisons: {nseq * (nseq - 1) // 2}")
print(f"Nucleotide diversity (pi): {pi}")
