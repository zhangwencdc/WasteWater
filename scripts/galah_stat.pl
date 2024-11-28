#!/usr/bin/perl
use strict;
use warnings;

#Usage:galah clusterÍ³¼Æ

my $file=$ARGV[0];#cluster.tsv
my %f;
open(F,$file);
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	$f{$a[0]}++;
}
close F;
open(O,">$file.stat");
my @k=sort {$f{$b}<=>$f{$a}} keys %f;
foreach my $k (@k) {
	print O "$k\t$f{$k}\n";
}