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

  foreach my $seq ($aln->each_seq()) {
    print $seq->display_id()." not multiple of 3! (".$seq->length().")\n" if ($seq->length() % 3 != 0);
  }

  #eval {
    my $stats = Bio::Align::DNAStatistics->new();
    my $result = $stats->calc_all_KaKs_pairs($aln);
    my ($Da, $Ds, $Dn, $N, $S, $S_d, $N_d);
    for my $an (@$result) {
      for (sort keys %$an ) {
        next if /Seq/;
        if($_ eq "D_n"){$Dn = $an->{$_}};
        if($_ eq "D_s"){$Ds = $an->{$_}};
        if($_ eq "S_d"){$S_d = $an->{$_};}
        if($_ eq "N_d"){$N_d = $an->{$_};}
        if($_ eq "S"){$S = $an->{$_};}
        if($_ eq "N"){$N = $an->{$_};}
      }
    }
    print Dumper($result);
    #print "Ka: $Dn\nKs: $Ds\n";
  #}
}
