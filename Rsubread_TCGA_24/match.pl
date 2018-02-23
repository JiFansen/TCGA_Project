#!/usr/bin/perl -w
use strict;
use warnings;
my %ref;
my @str;
my @sample;
my ($i,$j);
my $key;
while(<>)
{
	chomp;
	@str = split /\t/,$_;
	$ref{$str[0]} = $str[1];
}
open FILE,"./sample";
chomp(@sample = <FILE>);
close FILE;
foreach $i(0..$#sample)
{
	foreach $key(keys(%ref))
	{
		if($sample[$i] eq $key)
		{
			print "$sample[$i]\t$ref{$key}\n";
		}
	}
}
