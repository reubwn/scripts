#!/usr/bin/env perl

use strict;
use warnings;
use Sort::Naturally;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $usage = "
SYNOPSIS:
  Fetches and prints all downstream connections for an initial input scaffold name.
  I.e., also finds any connections to scaffolds connected to the input scaffold, so the whole connection path is returned.

OPTIONS:
  -c|--collinearity [FILE]   : collinearity file from MCScanX
  -s|--search       [STRING] : chromosome / scaffold name to get connections for
  -h|--help                  : print this message

USAGE:
  >> get_all_connections_from_collinarity.pl -c xyz.collinearity -s Ar1
\n";

my ($collinearity, $search, $help);

GetOptions (
  'collinearity|c=s' => \$collinearity,
  'search|s=s'       => \$search,
  'help|h'           => \$help,
);

die $usage if $help;
die $usage unless ($collinearity && $search);

my %connections;
open (my $IN, $collinearity) or die $!;

while (<$IN>) {
  if ($_ =~ m/(\w+\d+)\&(\w+\d+)\s/) {
    push @{ $connections{$1} }, $2;
    push @{ $connections{$2} }, $1;
  }
}

chomp $search;
my @all = @{$connections{$search}};
foreach ( nsort @{$connections{$search}} ) {
  push @all, @{$connections{$_}};
}

my @result = uniq(@all);
print join (",", (nsort @result)), "\n";

__END__
