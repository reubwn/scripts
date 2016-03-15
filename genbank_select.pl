#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
NAME
	genbank_select.pl

SYNOPSIS
	Reads in FEATURE annotations in Genbank files, outputs sequence data in fasta format.
	
	Features
	  - Always returns a fasta sequence, of either source (i.e., genomic) or gene/CDS
	  - Seq name is taken from the result of -t [TAGNAME] flag
	  - Note that some CDS don't have certain tags (e.g. gene names) and therefore would not get output if '-t gene' was specified
	  - To output protein translations, don't do '-t translation'; rather do '-o prot' (or '-o both')

OPTIONS
	-g|genbank : genbank format file (required)
	-p|ptag    : primary tag to read from (source, CDS) [default: CDS]
	-t|tag     : secondary tag to source within primary tag (e.g., gene, locus_tag, product, etc...) [default: locus_tag]
	-o|output  : output sequence type (nucl, prot, both) [default: nucl]
	--regexp   : perl regular expression - acts in conjunction with -t [TAGNAME], e.g. '-t gene --regexp \"xyl[A|B]\" -o prot' would pull out xylA and xylB translations from genes xylA and xylB
	-w|ins     : turn off case sensitivity for regexp flag, e.g., --regexp \"porin\" pulls out porin and Porin etc
	-i|include : file of seqnames to include, one per line. Assumes seqname is 1st column if there is other text in file. Note seqname MUST MATCH with -t [TAGNAME] provided or will return null
	-e|exclude : file of seqnames to exclude, same rules as above
	-d|delim   : delimiter for include|exclude file [default: any whitespace]
	-s|sort    : how to sort output seqs (natural, length, revlength, alpha, revalpha, random) [default: natural]
	-l|length  : length filter [default: none]

USAGE
	genbank_select.pl -g [GENBANKFILE] [options] > [FASTAFILE]

EXAMPLES
	genbank_select.pl -g genome.gbk	-p source -l 200 ## get all genomic sequences >=200
	genbank_select.pl -g genome.gbk -t gene --rexexp \"rpo[A-D]\" -o prot -s revlength ## get translations of rpo genes sorted smallest -> longest
	genbank_select.pl -g genome.gbk -p product --regexp \"amylase\" -w ## get CDSs of all genes with 'amylase' in protein description
\n";

## declare variables
my $gb_file;
my $prefix;
my $regexp;
my $insensitive;
my $includefile;
my %include_seqnames;
my $excludefile;
my %exclude_seqnames;
my $length;
my %fetched_seqs;
my $null;
my $help;

## set some defaults
my $p_tag = "CDS";
my $tag = "locus_tag";
my $output = "nucl";
my $numfasta = "all";
my $sort = "natural";
my $delimiter = "\\s+";
my $case = "u";

GetOptions (
	"g|gb|genbank=s" => \$gb_file,
	"p|ptag:s" => \$p_tag,
	"t|tag:s" => \$tag,
	"o|output:s" => \$output,
	"s|sort:s" => \$sort,
	"n|numfasta:i" => \$numfasta,
	"l|length:i" => \$length,
 	"prefix:s" => \$prefix,
 	"regexp:s" => \$regexp,
	"w|ins|insensitive" => \$insensitive,
 	"i|inc|include|includefile:s" => \$includefile,
 	"e|ex|exc|exclude|excludefile:s" => \$excludefile,
	"d|delimiter:s" => \$delimiter,
 	"c|case:s" => \$case,
	"h|help" => \$help,
);

die $usage if $help;
die $usage unless $gb_file;

## read gb file
my $seqio_object = Bio::SeqIO->new(-file => $gb_file, -format => "genbank");

## get include seqs if exists
if ($includefile){
	open (my $FH, $includefile) or die "$!\n";
	while (<$FH>){
		chomp;
		## splits on $delimiter, default is whitespace
		my @a = split $delimiter;
		## assumes seq names are in 1st column
		$include_seqnames{$a[0]} = ();
	}
	close $FH;
}

## iterate over seqs
while(my $seq_object = $seqio_object->next_seq){
	if ($p_tag eq 'CDS'){ ## get CDS seq features
		for my $feat_object ($seq_object->get_SeqFeatures) {
			if ($feat_object->primary_tag eq 'CDS'){
				if ($feat_object->has_tag($tag)){ ## check it has $tag tag
	
					## get $tag name and sequence data for each CDS $feat_obj within $seq_obj
					my %gene_seqs = map { $_->get_tag_values($tag), $_->spliced_seq->seq } $feat_object;
					my %protein_seqs = map { $_->get_tag_values($tag), $_->get_tag_values('translation') } $feat_object;

					## select depending on input
					if ($includefile){
						foreach (keys %include_seqnames){
							$fetched_seqs{"$_ [cds]"} = $gene_seqs{$_} if $output =~ m/(nucl|both)/ and $gene_seqs{$_};
		                                	$fetched_seqs{"$_ [translation]"} = $protein_seqs{$_} if $output =~ m/(prot|both)/ and $protein_seqs{$_};
						}
					} elsif ($excludefile){
						## TODO
					
					} elsif ($regexp){
						## apply regexp on $tag provided
						foreach (keys %gene_seqs){
							if ($insensitive){
								## turn off case sensitivity
								if ($_ =~ m/$regexp/i){
									$fetched_seqs{"$_ [cds]"} = $gene_seqs{$_} if $output =~ m/(nucl|both)/ and $gene_seqs{$_};
									$fetched_seqs{"$_ [translation]"} = $protein_seqs{$_} if $output =~ m/(prot|both)/ and $protein_seqs{$_};
								}
							} else {
								if ($_ =~ m/$regexp/){
                                                                 	$fetched_seqs{"$_ [cds]"} = $gene_seqs{$_} if $output =~ m/(nucl|both)/ and $gene_seqs{$_};
                                                                	 $fetched_seqs{"$_ [translation]"} = $protein_seqs{$_} if $output =~ m/(prot|both)/ and $protein_seqs{$_};
                                                        	 }
							}
						}
					} else {
						## print all seqs of primary tag CDS
						foreach (keys %gene_seqs){
        	                                	$fetched_seqs{"$_ [cds]"} = $gene_seqs{$_} if $output =~ m/(nucl|both)/ and $gene_seqs{$_};
                                                	$fetched_seqs{"$_ [translation]"} = $protein_seqs{$_} if $output =~ m/(prot|both)/ and $protein_seqs{$_};
                                        	}
					}
				} else {
					## sequence doesn't have $tag provided
					$null++; 
				}
			}
		}
	} elsif ($p_tag eq 'source'){ ## fetch raw sequence, under 'origin'
		$fetched_seqs{$seq_object->display_id} = $seq_object->seq;
	}
}

## apply length filter is set
if ($length){
	my %seq_lengths = map { $_, length($fetched_seqs{$_}) } keys %fetched_seqs;
	foreach (keys %seq_lengths){
		if ($seq_lengths{$_} < $length){
			delete $fetched_seqs{$_};
		}
	}
}

## print fetched seqs
if ($sort eq "random"){ ## output with no sort, which is kinda random
	foreach (keys %fetched_seqs){
		print "\>$_\n$fetched_seqs{$_}\n";
	}
} elsif ($sort eq "length"){ ## sort longest seq to shortest
	my %seq_lengths = map { $_, length($fetched_seqs{$_}) } keys %fetched_seqs;
	foreach (sort {$seq_lengths{$b} <=> $seq_lengths{$a}} keys %seq_lengths){
                print "\>$_\n$fetched_seqs{$_}\n";
        }
} elsif ($sort eq "revlength"){ ## sort shortest to longest
        my %seq_lengths = map { $_, length($fetched_seqs{$_}) } keys %fetched_seqs;
        foreach (sort {$seq_lengths{$a} <=> $seq_lengths{$b}} keys %seq_lengths){
                print "\>$_\n$fetched_seqs{$_}\n";
        }
} elsif ($sort eq "alpha"){ ## alphanumeric sort
	foreach (sort {$a cmp $b} keys %fetched_seqs){
                print "\>$_\n$fetched_seqs{$_}\n";
        }
} elsif ($sort eq "revalpha"){ ## reverse alphanumeric sort
	foreach (sort {$b cmp $a} keys %fetched_seqs){
                print "\>$_\n$fetched_seqs{$_}\n";
        }
} else {
	## default is natural sort (nicest)
	foreach (nsort keys %fetched_seqs){
                print "\>$_\n$fetched_seqs{$_}\n";
        }
}

print STDERR "[WARN]: $null sequences did not have tag \'$tag\' and therefore could not printed\n" if $null;
print STDERR "\n[INFO]: Printed ".(scalar(keys %fetched_seqs))." total sequences.\n\n";

__END__
