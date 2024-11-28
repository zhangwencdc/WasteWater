#!/usr/bin/perl
use strict;
use warnings;



my $data=$ARGV[0]; #PathogenCore marker基因的合并文件 大于1000bp

open(F,$data);my %anno;
my %species;
while(1){
		my $line=<F>;
		unless($line){last;}
		chomp $line;
		unless(substr($line,0,1) eq ">"){next;}
		my @b=split".fna_",$line;
		my $species=pop @b;
		$species{$species}++;
		$anno{substr($line,1)}=$species;
}
close F;

open(OUT,">$data.stat");
my @species=sort keys %species;
foreach my $sp (@species) {
	print OUT "$sp,$species{$sp}\n";
}