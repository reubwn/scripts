#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use Getopt::Long;
use Pod::Usage;
use Text::CSV;

# -----------------------------
# Defaults
# -----------------------------
my $input_dir    = "aligned_fastas";
my $pairwise_out = "pairwise_distances.csv";
my $summary_out  = "pi_summary.csv";
my $help         = 0;
my $man          = 0;

GetOptions(
    "input|i=s"    => \$input_dir,
    "pairwise|p=s" => \$pairwise_out,
    "summary|s=s"  => \$summary_out,
    "help|h"       => \$help,
    "man"          => \$man
) or pod2usage(2);

pod2usage(
    -verbose => 2,
    -exitval => 0
) if $man;

pod2usage(
    -verbose => 1,
    -exitval => 0
) if $help;

-d $input_dir or die "Error: input directory '$input_dir' does not exist\n";

my %valid = map { $_ => 1 } qw(A C G T);

# -----------------------------
# Find FASTA files
# -----------------------------
my @files;
find(
    sub {
        return unless -f $_;
        return unless $_ =~ /\.(fa|fasta|fas|fna)$/i;
        push @files, $File::Find::name;
    },
    $input_dir
);

die "No FASTA files found in '$input_dir'\n" unless @files;

# -----------------------------
# Open outputs
# -----------------------------
my $csv_pairwise = Text::CSV->new({ binary => 1, eol => "\n" })
  or die "Cannot create CSV object for pairwise output\n";

my $csv_summary = Text::CSV->new({ binary => 1, eol => "\n" })
  or die "Cannot create CSV object for summary output\n";

open my $pw_fh, ">", $pairwise_out
  or die "Cannot open '$pairwise_out' for writing: $!";
open my $sw_fh, ">", $summary_out
  or die "Cannot open '$summary_out' for writing: $!";

$csv_pairwise->print($pw_fh, [
    qw(file seq1 seq2 comparable_sites differences distance)
]);

$csv_summary->print($sw_fh, [
    qw(file n_sequences alignment_length n_pairs pi)
]);

# -----------------------------
# Process each file
# -----------------------------
for my $file (@files) {
    my ($ids_ref, $seqs_ref) = read_fasta_alignment($file);
    my @ids  = @$ids_ref;
    my @seqs = @$seqs_ref;

    my $n = scalar @ids;
    my $L = length($seqs[0]);

    my $sum_dist = 0;
    my $n_dist   = 0;
    my $n_pairs  = 0;

    for (my $i = 0; $i < $n - 1; $i++) {
        for (my $j = $i + 1; $j < $n; $j++) {
            my ($comparable, $differences, $distance) =
                pairwise_distance($seqs[$i], $seqs[$j], \%valid);

            $n_pairs++;

            $csv_pairwise->print($pw_fh, [
                $file,
                $ids[$i],
                $ids[$j],
                $comparable,
                $differences,
                defined $distance ? $distance : "NA"
            ]);

            if (defined $distance) {
                $sum_dist += $distance;
                $n_dist++;
            }
        }
    }

    my $pi = $n_dist ? ($sum_dist / $n_dist) : "NA";

    $csv_summary->print($sw_fh, [
        $file,
        $n,
        $L,
        $n_pairs,
        $pi
    ]);
}

close $pw_fh;
close $sw_fh;

# -----------------------------
# Subroutines
# -----------------------------

sub read_fasta_alignment {
    my ($file) = @_;

    open my $fh, "<", $file or die "Cannot open '$file': $!";

    my @ids;
    my @seqs;
    my ($current_id, $current_seq) = ("", "");

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/;

        if ($line =~ /^>(.*)/) {
            if ($current_id ne "") {
                push @ids,  $current_id;
                push @seqs, uc($current_seq);
            }
            $current_id  = $1;
            $current_id  =~ s/^\s+|\s+$//g;
            $current_seq = "";
        } else {
            $line =~ s/\s+//g;
            $current_seq .= $line;
        }
    }

    if ($current_id ne "") {
        push @ids,  $current_id;
        push @seqs, uc($current_seq);
    }

    close $fh;

    die "No sequences found in '$file'\n" unless @seqs;

    my $len = length($seqs[0]);
    for my $seq (@seqs) {
        die "Sequences are not aligned in '$file'\n"
            unless length($seq) == $len;
    }

    return (\@ids, \@seqs);
}

sub pairwise_distance {
    my ($seq1, $seq2, $valid_ref) = @_;

    my $len = length($seq1);
    my $comparable  = 0;
    my $differences = 0;

    for (my $k = 0; $k < $len; $k++) {
        my $b1 = substr($seq1, $k, 1);
        my $b2 = substr($seq2, $k, 1);

        next unless exists $valid_ref->{$b1} && exists $valid_ref->{$b2};

        $comparable++;
        $differences++ if $b1 ne $b2;
    }

    my $distance;
    if ($comparable > 0) {
        $distance = $differences / $comparable;
    } else {
        $distance = undef;
    }

    return ($comparable, $differences, $distance);
}

__END__

=head1 NAME

pairwise_pi.pl - Compute pairwise nucleotide distances and per-file pi from aligned FASTA files

=head1 SYNOPSIS

perl pairwise_pi.pl [options]

  --input,    -i   Input directory containing aligned FASTA files
  --pairwise, -p   Output CSV file for pairwise distances
  --summary,  -s   Output CSV file for per-file pi summary
  --help,     -h   Show this help message
  --man            Show full documentation (DESCRIPTION, EXAMPLES, etc.)

=head1 EXAMPLES

perl pairwise_pi.pl

perl pairwise_pi.pl 
  --input my_alignments
  --pairwise results_pairs.csv
  --summary results_summary.csv

=head1 DESCRIPTION

This script scans an input directory for FASTA files with extensions:
.fa, .fasta, .fas, .fna

For each aligned FASTA file, it computes pairwise nucleotide distances
using only sites where both sequences have A, C, G, or T.

Outputs:
- Pairwise results for each sequence pair
- Mean pairwise nucleotide diversity (pi) for each file

=cut
