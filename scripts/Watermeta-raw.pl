#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use Data::Dumper;
use File::Path;  
use Cwd;
my $path = getcwd;
#Usage:基于raw data的污水样本检测
##
my ($Keyname,$Ref,$input,$database,$type,$fa,$r,$level,$Outdir,$verbose,$db,$snp,$fq1,$fq2,$marker);
my ($Verbose,$Help);


GetOptions(
        "R:s"=>\$Ref, #参考基因组
        "L|level:s"=>\$level,  #默认Level: S  ,S代表Species，也可以设置为G，代表Genus
        "Out|O:s"=>\$Outdir,
        "verbose"=>\$Verbose,
        "Tag|Key:s" =>\$Keyname,
		"r|rlen:s" =>\$r,
		"database|DB:s"=>\$db,###比对数据库
		"1|Fq1:s"=>\$fq1,###必选 fq1
		"2|Fq2:s"=>\$fq2,###fq2
		"M|marker:s"=>\$marker,###marker
        "help"=>\$Help
);
####
$Keyname ||= "Input";
$Outdir ||= ".";
$type ||=1;
#$level ||="S";
$r ||=150;
$marker ||= $Bin."/Marker-v2.fasta";

die `pod2text $0` if ($Help);

if(substr($fq1,length($fq1)-3,3)=~/.gz/){system "gunzip $fq1\n";$fq1=substr($fq1,0,length($fq1)-3);}
if(defined $fq2){if(substr($fq2,length($fq2)-3,3)=~/.gz/){system "gunzip $fq2\n";$fq2=substr($fq2,0,length($fq2)-3);}}
###FastV
print "Fastv searching\n";
if(defined $fq2){
	system "fastv -i $fq1 -I $fq2 -o $Outdir/$Keyname.fq1.fastv.fq -O $Outdir/$Keyname.fq2.fastv.fq -c /home/zhangwen/Data/Fastv/microbial.kc.fasta -h $Outdir/$Keyname.fastv.html -j $Outdir/$Keyname.fastv.json\n";
}else{
	system "fastv -i $fq1 -o $Outdir/$Keyname.fq1.fastv.fq -c /home/zhangwen/Data/Fastv/microbial.kc.fasta -h $Outdir/$Keyname.fastv.html -j $Outdir/$Keyname.fastv.json\n";
}
system "rm -rf $Keyname.*.fastv.fq\n";
###Kraken2
print "Kraken2 searching\n";
if(defined $fq2){
	print "kraken2 --db /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE/ --report $Outdir/$Keyname.report --output $Outdir/$Keyname.kraken --paired $fq1 $fq2\n";
	system "kraken2 --db /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE/ --report $Outdir/$Keyname.report --output $Outdir/$Keyname.kraken --paired $fq1 $fq2\n";
	system "bracken -d /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE -i $Outdir/$Keyname.report -o $Outdir/$Keyname.report.bracken.species  -l S -r $r\n";
	system "bracken -d /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE -i $Outdir/$Keyname.report -o $Outdir/$Keyname.report.bracken.genus  -l G -r $r\n";
	system "perl /home/zhangwen/bin/kraken2-translate.pl $Outdir/$Keyname.report >$Outdir/$Keyname.report.txt\n";
	system "ktImportText $Outdir/$Keyname.report.txt -o $Outdir/$Keyname.kraken.html";
}else{
	system "kraken2 --db /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE/ --report $Outdir/$Keyname.report --output $Outdir/$Keyname.kraken $fq1 \n";
	system "bracken -d /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE -i $Outdir/$Keyname.report -o $Outdir/$Keyname.report.bracken.species  -l S -r $r\n";
	system "bracken -d /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE -i $Outdir/$Keyname.report -o $Outdir/$Keyname.report.bracken.genus  -l G -r $r\n";
	system "perl /home/zhangwen/bin/kraken2-translate.pl $Outdir/$Keyname.report >$Outdir/$Keyname.report.txt\n";
	system "ktImportText $Outdir/$Keyname.report.txt -o $Outdir/$Keyname.kraken.html";
}
system "rm $Outdir/$Keyname.kraken\n";
###VFF
print "Virulent factor searching\n";
if(defined $fq2){
	#system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -2 $fq2 -T /home/zhangwen/Data/VFF/VFDB_setA_nt.fas -o $Keyname.vff.setA\n ";
	system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -2 $fq2 -T /home/zhangwen/Data/VFF/VFDB_setB_nt.fas -o $Keyname.vff.setB\n ";
}else{
	#system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -T /home/zhangwen/Data/VFF/VFDB_setA_nt.fas -o $Keyname.vff.setA\n ";
	system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -T /home/zhangwen/Data/VFF/VFDB_setB_nt.fas -o $Keyname.vff.setB\n ";
}
system "rm $Keyname.*.sam* $Keyname.*.bam*\n";
##Resistance
print "Resistance genes searching\n";
if(defined $fq2){
	system "perl  /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -T /home/zhangwen/bin/Resistance/res.fas -1 $fq1 -2 $fq2 -o $Keyname.res \n ";
}else{
	system "perl  /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -T /home/zhangwen/bin/Resistance/res.fas -1 $fq1 -o $Keyname.res \n ";
}
#system "cp $Keyname.bam.sort.coverage $Keyname.res.report\n";
system "rm  $Keyname.*.sam* $Keyname.*.bam*\n";
###宿主鉴定
print "Host searching\n";
if(defined $fq2){
	system "perl  /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl  -1 $fq1 -2 $fq2 -o $Keyname.host -T /home/zhangwen/Data/Mitochondrion/Mitochondrion.fasta\n ";
}else{
	system "perl  /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl  -1 $fq1 -o $Keyname.host -T /home/zhangwen/Data/Mitochondrion/Mitochondrion.fasta\n ";
}
system "rm $Keyname.*.sam* $Keyname.*.bam*\n";
###Ref
if(defined $Ref){
	print "Target searching\n";
	if(defined $fq2){
	system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -2 $fq2 -T $Ref -o $Keyname.target\n ";
	}else{
		system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -T $Ref -o $Keyname.target\n ";
	}
system "rm $Keyname.*.sam* $Keyname.*.bam*\n";
}

##Marker 鉴定
print "Host searching\n";
if(defined $fq2){
	#system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -2 $fq2 -T /home/zhangwen/Data/VFF/VFDB_setA_nt.fas -o $Keyname.vff.setA\n ";
	system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -2 $fq2 -T $marker -o $Keyname.marker\n ";
}else{
	#system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -T /home/zhangwen/Data/VFF/VFDB_setA_nt.fas -o $Keyname.vff.setA\n ";
	system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $fq1 -T $marker -o $Keyname.marker\n ";
}
system "rm $Keyname.*.sam* $Keyname.*.bam*\n";
#system "multiqc $file -c /home/zhangwen/bin/multiqc_config.yaml\n";
print "Success Finished\n";
###############################
sub norRNA{
	my $fq=shift @_;
	my $kraken=shift @_;
	open(FILE,$kraken);
	my %rrna;
	while(1){
		my $line=<FILE>;
		unless ($line) {last;
		}
		chomp $line;
		my @a=split"\t",$line;
		if($a[0]eq "C"  ){$rrna{$a[1]}++;}
	}
	close FILE;
	my $result=$fq.".norRNA";
	open(OUT,">$result");
	open(FQ,$fq);
	while(1){
			my $line=<FQ>;
			unless($line){last;}
			chomp $line;
			my $seq=<FQ>;
			my $a=<FQ>;
			my $b=<FQ>;
			
			my @name=split" ",$line;
			my $name=substr($name[0],1);
			if(exists $rrna{$name}){next;}
			
			print OUT "$line\n$seq$a$b";
	}
	close FQ;
	return $result;
}

sub nohuman{
	my $fq=shift @_;
	my $kraken=shift @_;
	open(FILE,$kraken);
	my %human;
	while(1){
		my $line=<FILE>;
		unless ($line) {last;
		}
		chomp $line;
		my @a=split"\t",$line;
		if($a[2]eq "9605" || $a[2]eq "9606" ){$human{$a[1]}++;}
	}
	close FILE;
	my $result=$fq.".nohuman";
	open(OUT,">$result");
	open(FQ,$fq);
	while(1){
			my $line=<FQ>;
			unless($line){last;}
			chomp $line;
			my $seq=<FQ>;
			my $a=<FQ>;
			my $b=<FQ>;
			
			my @name=split" ",$line;
			my $name=substr($name[0],1);
			if(exists $human{$name}){next;}
			
			print OUT "$line\n$seq$a$b";
	}
	close FQ;
	return $result;
}