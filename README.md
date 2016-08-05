# Scripts
Useful scripts for various things.

## n50.pl
Calculate various assembly-like metrics for fasta files, including span, number of sequences, n50, n90, %NNNs, %GC... etc.

### Usage
```perl
n50.pl genome1.fasta [genome2.fasta...]
```

## genbank_select.pl
Extract and filter sequences from a genbank-format file. Sequences can be genomic (source), CDS or translations (protein). In development.

### Usage
Type `genbank_select -h` for options.

## bam2tab.pl
Script for filtering matepair reads from a SAM/BAM file.

It optionally removes any pair of reads that (a) map to the same contig (since these have no scaffolding information), and (b) pairs where either read maps to within a specified distance from the edge of a contig. This is because there can be a lot of contamination from short-insert, FR reads in an MP library, which may lead to mis-scaffolds. By removing any read mapping to within X bases of a contig edge, you can ensure that pairs which pass the filter have a _minimum mapping distance_ which is some value greater than the insert size of the contaminating FR reads.

Edge distance is calculated as:

```
contig2
==========+++++++++++++
  3'<----5'
    read2(R)            contig1
                       ++++++++==================
                              5'---->3'
                                read1(F)
the regions marked +++ must be > $edge_distance
```
when read1 is on the + strand and

```
contig1
==========+++++++++++++
  3'<----5'
    read1(F)            contig2
                        ++++++++==================
                              5'---->3'
                                   read2(R)
the regions marked +++ must be > $edge_distance
```
when read1 is on the - strand.

### Usage
Type `bam2tab.pl -h` for options:

```
OPTIONS:
  -i|--in            [FILE]  : SAM/BAM file [required]
  -f|--fasta         [FILE]  : fasta file of contigs [required]
  -d|--edge_distance [INT]   : filter based on minimum mapping distance from contig edge [default: 0]
  -s|--same_contig           : filter reads mapping to the same contig (i.e., are redundant) [default: no]
  -o|--output        [FILE]  : output file name [default: mapping_table.txt]
  -t|--outtype       [STR]   : Output format: 'table' => table [default], 'sam' => sam format
  -h|--help                  : prints this help message
```

### Outputs
Either TAB formatted for use in SSPACE or as SAM.
