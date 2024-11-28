#!/usr/bin/perl
use strict;
#use warnings;
use File::Basename qw(basename dirname);
#/home/zhangwen/project/2024Water/Analysis/Bacteria/Bracken/*.report.bracken.genus

my $sample=shift @ARGV;#../Sample.info
my @file=@ARGV;
my $cut=0.0001;#0.0001 仅考虑丰度超过万分之一的菌


#读取sample.info
open(F,$sample)||die;

my $l=<F>;
chomp $l;my %anno;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	$anno{$a[1]}=$l;# 需根据sample 顺序调整
	#print "$sample,$a[1]\n";
}
close F;
my %target;

foreach my $f (@file) {
	open(F,$f);
	my $n=basename($f);
	my @a=split"\\.",$n;
	my $name=$a[0];
#print "$f,$name,$anno{$name}\n";
	unless(exists $anno{$name}){next;}
	#	print "$f,$name\n";
	while(1){
		my $l=<F>;
		unless($l){last;}
		chomp $l;
		my @a=split"\t",$l;
	
		unless($a[6]>=$cut){next;}
		$target{$a[0]}{$name}=$a[6];
		#print "$l\n";
	}
	close F;
}

my @target=sort keys %target;
print "Target,";
foreach my $t (@target) {
	print "$t,";
}
print "Anno\n";
my @sample=sort keys %anno;
foreach my $sample (@sample) {
	print "$sample,";
	foreach my $t (@target) {
		print "$target{$t}{$sample},";
	}
	print "$anno{$sample}\n";
}