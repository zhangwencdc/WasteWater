#!/usr/bin/perl
#use strict;
#use warnings;

#Usage:Report

my $sample=$ARGV[0];#sample.info 样本信息表
my $target=$ARGV[1];#target.list 拟调查的genus列表
my $stat=$ARGV[2];#Euk.cat.50 bowtiereport合并结果

open(F,$sample)||die;
my %genome;my %country;my %year;my %country_n;my %year_n;
my $l=<F>;
chomp $l;my %sample;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	$genome{$a[1]}=$a[6];
	$country{$a[1]}=$a[7];
	$year{$a[1]}=$a[8];
	$sample{"Sum"}++;
	$sample{$a[7]}{"Country"}++;
	$sample{$a[8]}{"Year"}++;
	$country_n{$a[7]}++;
	$year_n{$a[8]}++;
}
close F;
my @id=sort keys %sample;
open(FILE,$stat)||die;
my $l=<FILE>;chomp $l;my %pos;my %target; my %pos_country;my %pos_year;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
		my @c=split">",$l;
		my $anno=pop @c;
		my @d=split" ",$anno;
		my $genus=$d[1];
		my @a=split"\\.",$l;
	#	print "$genus\n";
			unless(exists $genome{$a[0]}){next;}
		$pos{$a[0]}{$genus}++;

		print "$l\n";
		$target{$genus}++;
		
		my $id=$country{$a[0]};
			if(exists $sample{$id}{"Country"}){
				$pos_country{$id}{$genus}++;
			}
			$id=$year{$a[0]};
			if(exists $sample{$id}{"Year"}){
				$pos_year{$id}{$genus}++;
			}
	
}
close FILE;

open(T,$target);
my @sample_num=keys %genome; my $sample_num=@sample_num; my %t;
print "Target,Positive sample,Total sample,Percentage\n";
while(1){
	my $l=<T>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	#print "$l,$target{$l}\n";
	if(exists $target{$l}){my $p=$target{$l}/$sample_num;
	print "$l,$target{$l},$sample_num,$p\n";$t{$l}++;
	}
}
close T;
#die;
my @country=sort keys %country_n;
print "Country,"; 
my @target=sort keys %t;
foreach my $s (@country) {
	print "$s,";
}
print "\nNum,";
foreach my $s (@country) {
	print "$country_n{$s},";
}
print "\n";
foreach my $target (@target) {
	print "$target,";
	foreach my $s (@country) {
		print "$pos_country{$s}{$target},";
	}
	print "\n";
}

my @year=sort keys %year_n;
print "Year,"; 
my @target=sort keys %t;
foreach my $s (@year) {
	print "$s,";
}
print "\nNum,";
foreach my $s (@year) {
	print "$year_n{$s},";
}
print "\n";
foreach my $target (@target) {
	print "$target,";
	foreach my $s (@year) {
		print "$pos_year{$s}{$target},";
	}
	print "\n";
}