#!/usr/bin/perl
#use strict;
#use warnings;
use File::Basename qw(basename dirname);
#Usage:Marker Report

my $sample=$ARGV[0];#sample.info 样本信息表
my $stat=$ARGV[1];#target.list Bacteria_marker.stat


open(F,$sample)||die;
my %genome;my %country;my %year;my %country_n;my %year_n;my %fasta;
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
	$fasta{$a[6]}=$a[1];
}
close F;
open(FILE,$stat)||die;
my $l=<FILE>;chomp $l;my %target;my %species;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my $name=basename($a[0]);
	unless(exists $fasta{$name}){next;}
	unless($a[4]>=0.05){next;}
	my @b=split"_",$a[1];
	my $target=$b[0]."_".$b[1];
	$target{$name}{$target}++;
	$species{$target}++;
	#print "$name,$target\n";
}
close FILE;

my @sample=sort keys %genome;
my @species=sort keys %species;
print "Sample,Name,Country,Year,"; my %c_sp;my %y_sp;
foreach my $sp (@species) {
	print "$sp,";
}
print "\n";
foreach my $sample (@sample) {
	my $sname=$genome{$sample};
	print "$sample,$sname,$country{$sample},$year{$sample},";
	foreach my $sp (@species) {
		if(exists $target{$sname}{$sp}){print "1,";
			$c_sp{$sp}{$country{$sample}}++;
			$y_sp{$sp}{$year{$sample}}++;
		}else{print ",";}
	}
	print "\n";
}

my @country=sort keys %country_n;
print "Country,";
foreach my $country (@country) {
	print "$country,";
}
print "\n";
print "Num,";
foreach my $country (@country) {
	print "$country_n{$country},";
}
print "\n";
foreach my $sp (@species) {
	print "$sp,";
	foreach my $country (@country) {
		print "$c_sp{$sp}{$country},";
	}
	print "\n";
}

my @year=sort keys %year_n;
print "year,";
foreach my $year (@year) {
	print "$year,";
}
print "\n";
print "Num,";
foreach my $year (@year) {
	print "$year_n{$year},";
}
print "\n";
foreach my $sp (@species) {
	print "$sp,";
	foreach my $year (@year) {
		print "$y_sp{$sp}{$year},";
	}
	print "\n";
}