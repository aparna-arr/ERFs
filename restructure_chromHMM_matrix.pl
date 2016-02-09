#!/usr/bin/env perl
use warnings;
use strict;

my ($input, $output) = @ARGV;
die "usage: $0 <chromHMMfile> <outfile>\n" unless @ARGV;

open (IN, "<", $input) or die "could not open $input\n";

my %states;
my @rows;
while(<IN>)
{
	my $line = $_;
	chomp $line;
	
	next if ($line =~ /^track/ || $line =~ /^\#/);

	my ($chr, $start, $end, $state, $trash) = split(/\t/, $line);

	$states{$state}{avg} += $end - $start;
	$states{$state}{count}++;
}

close IN; 

open (OUT, ">", $output);
my $string1;
my $string2;

foreach my $state (keys %states)
{
	$string1 .= "$state\t";
	$string2 .= ($states{$state}{avg} / $states{$state}{count}) . "\t";		
}

print OUT $string1 . "\n";
print OUT $string2 . "\n";

close OUT;

