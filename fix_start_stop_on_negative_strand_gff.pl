#!/usr/bin/env perl

## author: reubwn Aug 2020

use strict;
use warnings;

use Getopt::Long;

my $usage = "
SYNOPSIS:
  Somehow start_codon and stop_codon features are not reported correctly after GAG.
  Fix it.

OPTIONS:
  -g|--gff   [FILE] : Input GFF file [required]
  -o|--out   [FILE] : Output filename [fixed.gff]
  -h|--help         : Prints this help message
\n";

my ($gff_file,$help,$debug);
my $gff_outfile = "fixed_strand.gff";

GetOptions (
  'g|gff=s'     => \$gff_file,
  'o|outfile:s' => \$gff_outfile,
  'h|help'      => \$help,
  'd|debug'     => \$debug
);

die $usage if $help;
die $usage unless ($gff_file);

## IO
my %strand_hash;
open (my $IN1, $gff_file) or die $!;

## recurse GFF once
print STDERR "[INFO] Building strand hash\n";
while (my $line = <$IN1>) {
  chomp $line;
  if ($line =~ m/^#/) {
    next;
  }
  $line =~ s/;$//; ## trim trailing ';' if present
  my @F = split (/\s+/, $line);

  ## hash of mRNA IDs and their strand value
  if ( $F[2] eq "mRNA" ) {
    my $ID = $1 if ($F[8] =~ m/ID=(.+?)(;|$)/); ## inherit ID from ID
    $strand_hash{$ID} = $F[6];
  }
}
close $IN1;

## IO
open (my $IN2, $gff_file) or die $!;
open (my $OUT, ">$gff_outfile") or die $!;
print STDERR "[INFO] Writing to '$gff_outfile'\n";

## recurse GFF second time to make the edits
print STDERR "[INFO] Checking GFF...\n";
while (my $line = <$IN2>) {
  chomp $line;
  if ($line =~ m/^#/) {
    print $OUT "$line\n";
    next;
  }
  $line =~ s/;$//; ## trim trailing ';' if present
  my @F = split (/\s+/, $line);

  ## check start_codon
  if ( $F[2] eq "start_codon" ) {
    my $ID = $1 if ($F[8] =~ m/Parent=(.+?)(;|$)/); ## inherit ID from Parent
    if ($F[6] ne $strand_hash{$ID}) {
      print STDERR "[INFO] Conflict found: strand for $ID:start is '$F[6]' but should be '$strand_hash{$ID}'\n" if ( $debug );
      print $OUT join ("\t", @F[0..5], $strand_hash{$ID}, @F[7..$#F]) . "\n";
    } else {
      print $OUT "$line\n";
    }
  } elsif ( $F[2] eq "stop_codon" ) {
    my $ID = $1 if ($F[8] =~ m/Parent=(.+?)(;|$)/); ## inherit ID from Parent
    if ($F[6] ne $strand_hash{$ID}) {
      print STDERR "[INFO] Conflict found: strand for $ID:stop is '$F[6]' but should be '$strand_hash{$ID}'\n" if ( $debug );
      print $OUT join ("\t", @F[0..5], $strand_hash{$ID}, @F[7..$#F]) . "\n";
    } else {
      print $OUT "$line\n";
    }
  } else {
    print $OUT "$line\n";
  }
}
close $IN2;
close $OUT;

print STDERR "[INFO] Done " . `date`;

__END__
