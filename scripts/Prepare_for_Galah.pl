#!/usr/bin/perl
use strict;
use warnings;

#Usaage:prepare for Galah

my $checkv=$ARGV[0];#CheckV report  CheckV_out/quality_summary.tsv

#my $genomad=$ARGV[1];#Aquatic_wastewater_genomad_virus_summary.tsv

my $virus=$ARGV[1];#Aquatic_wastewater_genomad_virus.fna

my $outdir=$ARGV[2];#/home/zhangwen/project/2024Water/Analysis/Noval_Virus/Virus_seq

my %quality;
open(F,$checkv);
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless($a[1]>=1500){next;}  ##保留长度超过1500bp的序列
	if($a[7]=~/Low-quality/ || $a[7]=~/Not-determined/ || $a[7]=~/checkv_quality/){next;}
	print "$outdir/$a[0].fasta\n";
	system "seqkit grep -n -p $a[0] $virus>$outdir/$a[0].fasta\n";
}
close F;