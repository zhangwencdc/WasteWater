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

my @file=glob "$ARGV[1]/*.fasta";#/data/zhangwen/2023Aquatic/Data/Assemble
	open(OUT,">Marker.blat");
	print OUT "Filename\tMarker\tSeqID\tAlignLength\tGeneLength\n";
foreach my $file (@file) {
	system "blat $data $file tmp.blat";
	open(F,"tmp.blat");
	
	while(1){
		my $line=<F>;
		unless($line){last;}
		chomp $line;
		my @a=split"\t",$line; 
		unless($a[0]>0){next;}
		if($a[0]>=$a[14]*0.3 || $a[0]>=300){print OUT "$file\t$a[13]\t$a[9]\t$a[0]\t$a[14]\n";}
	}
	close F;
	system "rm -rf tmp.blat\n";
}

open(F,"Marker.blat");my %match;my %tmp;
while(1){
		my $line=<F>;
		unless($line){last;}
		chomp $line;
		my @a=split"\t",$line;
		
		unless($a[3]>=0.5*$a[4]){next;} #仅考虑align大于50%的比对结果
		if(exists $tmp{$a[0]}{$a[1]}){next;}
		$tmp{$a[0]}{$a[1]}++;
		my $anno=$anno{$a[1]};
		$match{$a[0]}{$anno}++;
}
close F;

open(O,">Marker.blat.stat");
my @target=sort keys %match;
my @species=sort keys %species;
foreach  my $target(@target) {
	foreach my $species (@species) {
		if(exists $match{$target}{$species}){
			my $per=$match{$target}{$species}/$species{$species};
			print O "$target\t$species\t$species{$species}\t$match{$target}{$species}\t$per\t";
			if($per>=0.05){print O "Positive\n";}else{print O "\n";}
		}
	}
}
close O;