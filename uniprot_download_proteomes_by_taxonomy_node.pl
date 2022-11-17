#!/usr/bin/env perl

use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Date;

# Taxonomy identifier of top node for query, e.g. 2 for Bacteria, 2157 for Archea, etc.
# (see https://www.uniprot.org/taxonomy)
my $top_node = $ARGV[0];

my $agent = LWP::UserAgent->new;

# Get a list of all reference proteomes of organisms below the given taxonomy node.
my $query_list = "https://rest.uniprot.org/proteomes/stream?query=reference:true+taxonomy_id:$top_node&format=list";

my $response_list = $agent->get($query_list);
die 'Failed, got ' . $response_list->status_line . ' for ' . $response_list->request->uri . "\n" unless $response_list->is_success;

# For each proteome, mirror its set of UniProt entries in compressed FASTA format.
for my $proteome (split(/\n/, $response_list->content)) {
  my $file = $proteome . '.fasta.gz';
  my $query_proteome = "https://rest.uniprot.org/uniprotkb/stream?query=proteome:$proteome&format=fasta&compressed=true";
  my $response_proteome = $agent->mirror($query_proteome, $file);

  if ($response_proteome->is_success) {
    my $release = $response_proteome->header('x-uniprot-release');
    my $date = $response_proteome->header('x-uniprot-release-date');
    print "File $file: downloaded entries of UniProt release $release ($date)\n";
  } elsif ($response_proteome->code == HTTP::Status::RC_NOT_MODIFIED) {
    print "File $file: up-to-date\n";
  } else {
    die 'Failed, got ' . $response_proteome->status_line . ' for ' . $response_proteome->request->uri . "\n";
  }
}
