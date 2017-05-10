#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq::EncodedSeq;
use Bio::Align::DNAStatistics;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Data::Dumper;

my @files = <./*$ARGV[0]>;
foreach (@files) {
  print "$_\n";
  my $read_aln = Bio::AlignIO -> new(-file=>$_, -format=>"fasta");
  my $aln = $read_aln->next_aln(); ##read aln

  my @seq_names;
  foreach my $seq ($aln->each_seq()) {
    if ($seq->length() % 3 != 0) {
      print $seq->display_id()." not multiple of 3! (".$seq->length().")\n";
    }
    push (@seq_names, $seq->display_id());
  }

  # ## check if all CDS seqs are multiple of 3:
  # foreach (keys %cds_seqs) {
  #   if ($cds_seqs{$_}->length() % 3 != 0) {
  #     print STDERR "\n[WARN] Seq ".$cds_seqs{$_}->display_id()." is not a multiple of 3; will trim ".($cds_seqs{$_}->length() % 3)." bases from 3' end\n";
  #     my $trimmed = $cds_seqs{$_}->subseq(1,(($cds_seqs{$_}->length()) - ($cds_seqs{$_}->length() % 3))); ##trims remainder off 3' end; returns a STRING, annoyingly
  #     $cds_seqs{$_} = Bio::Seq->new( -display_id => $_, -seq => $trimmed ); ##replace old seq with trimmed seq
  #     $trimmed_seqs++;
  #   }
  # }

  #eval {
    my $stats = Bio::Align::DNAStatistics->new();
    my @result = @{$stats->calc_KaKs_pair($aln, $seq_names[0], $seq_names[1])};
    # my %result = %{$r};
    if (exists($result[0]{D_n})) {
      print "Ka: $result[0]{D_n}\n";
    } else { print "Ka: -2\n"; }
    if (exists($result[0]{D_s})) {
      print "Ks: $result[0]{D_s}\n";
    } else { print "Ks: -2\n"; }

    # for my $an (@$result) {
    #   for (sort keys %$an ) {
    #     next if /Seq/;
    #     if($_ eq "D_n"){$Dn = $an->{$_}};
    #     if($_ eq "D_s"){$Ds = $an->{$_}};
    #     if($_ eq "S_d"){$S_d = $an->{$_};}
    #     if($_ eq "N_d"){$N_d = $an->{$_};}
    #     if($_ eq "S"){$S = $an->{$_};}
    #     if($_ eq "N"){$N = $an->{$_};}
    #   }
    # }
    print Dumper \@result;
    #print "Ka: $Dn\nKs: $Ds\n";
  #}
}
