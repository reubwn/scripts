#!/usr/bin/env perl

## author: reubwn May 2019

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;

my $usage = "
SYNOPSIS:
  Strips trailing asterisks (stop codons; '*') from protein sequences
  Also removes sequences with internal stop codons by default
  Can read from gzip

OPTIONS:
  -f|--fasta [FILE]   : fasta file [required]
  -s|--stop  [STRING] : set stop string (default '*')
  -i|--inplace        : do in-place replacement of file (with '*.bak' backup)
  -h|--help           : prints this help message
\n";

my ($fasta, $stop, $inplace, $help);

GetOptions (
  'f|fasta=s' => \$fasta,
  's|stop:s'  => \$stop,
  'i|inplace' => \$inplace,
  'h|help'    => \$help
);

die $usage if $help;
die $usage unless ($fasta);

my $search_string = quotemeta ($stop);

my (%stripped, %removed);

my $in;
if ($fasta =~ m/gz$/) { ## read from gzip
  $in = Bio::SeqIO -> new ( -file => "gunzip -c $fasta |", -format => "fasta" );
} else {
  $in = Bio::SeqIO -> new ( -file => $fasta, -format => "fasta" );
}

if ($inplace) { ## do in-place replacement of file
  print STDERR "[INFO] Doing in-place file replacement (backup saved to $fasta.bak)\n";
  ## save a backup
  `cp $fasta $fasta.bak`;

  open (my $TMP, ">$fasta.tmp") or die $!;
  while (my $seq_obj = $in->next_seq()) {
    my $seq_string = $seq_obj->seq();
    if ($seq_string =~ m/$search_string$/) { ## has terminal stop codon
      $seq_string =~ s/$search_string$//; ## remove it
      if ($seq_string =~ m/$search_string/) { ## also has internal stop codon
        $removed{$seq_obj->display_id()}++; ## dont print it
      } else {
        print $TMP ">" . $seq_obj->display_id() . "\n" . $seq_string . "\n";
        $stripped{$seq_obj->display_id()}++;
      }

    } elsif ($seq_string =~ m/$search_string/) { ## stop codon is internal; don't print
      $removed{$seq_obj->display_id()}++; ## dont print it
    } else {
      print $TMP ">" . $seq_obj->display_id() . "\n" . $seq_string . "\n";
    }
  }
  ## do the replacement
  `mv $fasta.tmp $fasta`;

} else { ## print to STDOUT
  while (my $seq_obj = $in->next_seq()) {
    my $seq_string = $seq_obj->seq();
    if ($seq_string =~ m/$search_string$/) { ## has terminal stop codon
      $seq_string =~ s/$search_string$//; ## remove it
      if ($seq_string =~ m/$search_string/) { ## also has internal stop codon
        $removed{$seq_obj->display_id()}++; ## dont print it
      } else {
        print STDOUT ">" . $seq_obj->display_id() . "\n" . $seq_string . "\n";
        $stripped{$seq_obj->display_id()}++;
      }

    } elsif ($seq_string =~ m/$search_string/) { ## stop codon is internal; don't print
      $removed{$seq_obj->display_id()}++; ## dont print it
    } else {
      print STDOUT ">" . $seq_obj->display_id() . "\n" . $seq_string . "\n";
    }
  }

}

print STDERR "[INFO] Number of seqs with trailing '$stop' stripped: " . (keys %stripped) . "\n";
print STDERR "[INFO] Number of seqs with internal '$stop' removed: " . (keys %removed) . "\n";
print STDERR "[INFO] Done " . `date`;
