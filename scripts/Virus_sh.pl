#!/usr/bin/perl
use strict;
use warnings;

my $file=$ARGV[0];#Sample.info
my $target="/home/zhangwen/Data/2023Genome_Database/DRMVirus/DRMVirus_2023.fasta";
my $script="/home/zhangwen/project/2024Water/Analysis/Virus/Target_bowtie-dell-Micro-v2.pl";
open(F,$file);
my $l=<F>;chomp $l;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my $fq1="/home/zhangwen/project/2024Water/Data/".$a[0]."/Raw/".$a[4];
	my $fq2="/home/zhangwen/project/2024Water/Data/".$a[0]."/Raw/".$a[5];
	print "$a[1]\n";
	if($fq1=~/[0-9a-zA-Z]/){
		if($fq2=~/[0-9a-zA-Z]/){
			system "perl $script -1 $fq1 -2 $fq2 -T $target -o $a[1]\n";
		}else{
			system "perl $script -1 $fq1  -T $target -o $a[1]\n";
		}
	}
}
close F;
