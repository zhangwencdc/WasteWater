#!/usr/bin/perl
use strict;
#use warnings;
use File::Basename qw(basename dirname);
#Usage: 合并bowtie report结果
#perl Bowtie_stat.pl /home/zhangwen/project/2024Water/Analysis/Euk 5 Euk.cat
my @file=glob "$ARGV[0]/*.report"; 
my $cutoff=$ARGV[1];# cutoff 5 指align 5%
open(O,">$ARGV[2]");
print O "Sample,Target,Target Length,Match Read Num,Align Length,Coverage%,Depth,Annotation\n";
my %genus;
foreach my $file (@file) {
	my $fname=basename($file);
	open(F,$file);
	while(1){
		my $l=<F>;
		unless ($l) {last;
		}
		chomp $l;
		my @a=split",",$l;
		my @b=split " ",$a[4];
		unless($b[0]>=$cutoff){next;}
		print O "$fname,$l\n";
		my @c=split">",$l;
		my $anno=pop @c;
		my @d=split" ",$anno;
		my $genus=$d[1];
		$genus{$d[1]}++;
	}
	close F;
}
close O;
my @genus=sort keys %genus;
print "Genus,Positive sample num\n";
foreach my  $g(@genus) {
	print "$g,$genus{$g}\n";
}