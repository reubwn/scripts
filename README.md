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

## bam_edgeFilter.pl
Script for filtering matepair reads from a SAM/BAM file.

It optionally removes any pair of reads that (a) map to the same contig (since these have no scaffolding information), and (b) pairs where either read maps to within a specified distance from the edge of a contig. This is because there can be a lot of contamination from short-insert, FR reads in an MP library that may lead to mis-scaffolding. By removing any read mapping to within X bases of a contig edge, you can ensure that pairs which pass the filter have a __minimum mapping distance__ which is some value greater than the upper bounds of the insert size distribution of the contaminating FR reads.

There are 4 edge distances for each mapped read:

```
      contig         read
start|++++++++++++<========++++++++++++++++++++++|end
      <----(1)---><------------(2)-------------->
      <--------(3)--------><--------(4)--------->
```

Thus 8 when both reads are considered; all 8 must be > threshold to pass filter.

### Usage
Type `bam_edgeFilter.pl -h` for options:

```
OPTIONS:
  -i|--in                [FILE] : SAM/BAM file [required]
  -f|--fasta             [FILE] : fasta file of contigs [required]
  -d|--min_edge_distance [INT]  : discard reads with edge distance < INT [default: 500]
  -s|--same_contig              : discard reads mapping to the same contig [default: no]
  -o|--out               [STR]  : output prefix [default: edgeFilter_output]
  -t|--outtype           [STR]  : output format: 'sam' => sam format [default], 'table' => table
  -h|--help                     : prints this help message
```

Recommend to set -d to __at least__ the average insert size of the contaminating FR reads.

### Outputs
Either as SAM (recommended) or as TAB formatted for use in SSPACE (untested).

### Pipeline
1. Map trimmed MP reads to assembly
2. Estimate insert size distributions of FR, RF, Tandem reads (e.g., with PicardTools CollectInsertSizeMetrics)
3. Run bam_edgeFilter.pl with -d <mean_FR_insert_size>
4. Convert output SAM file to fastq (PicardTools SamToFastq)
5. Scaffold contigs using filtered reads (e.g., SSPACE)

## bam2tab.pl
Old SAM filtering code. Obsoleted by bam_edgeFilter.pl.

## orthofinder_selectGroups.pl
Script for pulling out groups of interest from an orthofinder output file OrthologousGroups.txt, and writing each OG to a fasta file (eg., for downstream alignment etc).

Construct the 'select' string like so: \"ID1=X,ID2=Y,IDN=Z\", where ID1 is the unique species ID corresponding to the entity on the lefthandside of the delimiter given by -d (default is '|'), and X is the number of members from species 1 you want to select groups containing.

### Usage
```
  -i|--in     [FILE] : OrthologousGroups.txt file [required]
  -f|--fasta  [FILE] : path to directory of fasta files used to construct OrthologousGroups.txt [required]
  -d|--delim  [STR]  : delimiter used in the protein naming structure [assumes an OrthoMCL-style schema of "SPECIES_ID|PROTEIN_ID"]
  -s|--select [STR]  : select string, eg. "ID1=1,ID2=1,IDN=1" would return 1:1 orthologous groups
  -n|--noseqs        : Don't print sequences; just count the groups
  -o|--out    [STR]  : output prefix [default: selectGroups_output]
  -h|--help          : prints this help message
```

### Outputs
A file "selectGroups_output.groups.txt" with only selected groups written, and a directory "selectGroups_outputDir/" with fasta files containing sequence data of selected groups.

### Example
Select ONLY those groups where species ARIC and AVAG have 4 members each, and species RMAG and RMAC have 2 members each:
```
>> ./orthofinder_selectGroups.pl -i OrthologousGroups.txt -f fastas/ -s "ARIC=4,AVAG=4,RMAG=2,RMAC=2"
```
