#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#Usage:galah clusterÍ³¼Æ

my $file=$ARGV[0];#cluster.tsv
my $anno=$ARGV[1];#Aquatic_wastewater_genomad_virus_summary.tsv

open(FILE,$anno);
my %anno;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	$anno{$a[0]}=$a[10];
}

close FILE;


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
	my $name=basename($k);
	$name=substr($name,0,length($name)-6);
	print O "$k\t$f{$k}\t$anno{$name}\n";
}