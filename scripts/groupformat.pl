#!/usr/bin/perl
use strict;
#use warnings;
#≈‰∫œmeta.R π”√

my $group=$ARGV[0];#Group.txt

open(Group,$group);my %group;
my $l=<Group>;chomp $l;

while(1){
	my $l=<Group>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;	
	my $name=shift @a;
	foreach my $a (@a) {
		if($a=~/[0-9a-zA-Z]/ || $a eq "NA" || $a eq "na"){
		$group{$a}++;
		}
	}
	
}
close Group;
my @key=keys %group;
my $num=0;my %key;
foreach my $key (@key) {
	$key{$key}=$num;
	$num++;
}

open(Group,$group);
open(OUT,">group_format.txt");
my $l=<Group>;chomp $l;
my @a=split"\t",$l;	
foreach my $a (@a) {
	print OUT "$a,";
}
print OUT "\n";
while(1){
	my $l=<Group>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;	
	my $name=shift @a;
	print OUT "$name,";
	foreach my $a (@a) {
		if(exists $key{$a}){print OUT "$key{$a},";}else{print OUT ",";}
	}
	print OUT "\n";
}