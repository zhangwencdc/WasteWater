#!/usr/bin/perl
use strict;
#use warnings;
#≈‰∫œmeta.R π”√



my $matrix=$ARGV[0];#bray matrix;
my $group=$ARGV[1];
my $id=$ARGV[2];
print "Input: $matrix $group $id\n";

open(Group,$group);my %group;
my $l=<Group>;chomp $l;
while(1){
	my $l=<Group>;
	unless($l){last;}
	chomp $l;
	my @a=split",",$l;
	unless($a[$id]=~/[0-9a-zA-Z]/){next;}
	if($a[$id] eq "NA" || $a[$id] eq "na"){next;}
	$group{$a[0]}=$a[$id];
}
close Group;


my @group=keys %group;


open(FILE,$matrix);

my $line=<FILE>;chomp $line;
my @name=split",",$line;my $id=0;
open (OUT,">dis.csv");
open (O,">dis_v2.csv");
print OUT "ID,Sample1 name,Sample1 Group,Sample2 Name,Sample2 Group,Type,Bray\n";
print O "ID,Sample1 name,Sample1 Group,Sample2 Name,Sample2 Group,Type,Bray\n";
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	my @a=split",",$line;
	my $n=@a;
	foreach (1..($n-1)) {
		my $b=$name[$_];
		if($a[0] eq $b){next;}
		my $ga=$group{$a[0]};
		my $gb=$group{$b};
		my $dis=$a[$_];
		unless($ga=~/[0-9a-zA-Z]/ && $gb=~/[0-9a-zA-Z]/){next;}
		$id++;
		print OUT "$id,$a[0],$ga,$b,$gb,";my $type; my $type2;
		if($ga eq $gb){$type="Within $ga"; $type2="Within";}else{$type="Between";$type2="Between";}
		print OUT "$type,$dis\n";print O "$id,$a[0],$ga,$b,$gb,$type2,$dis\n";
	}
}
close FILE;