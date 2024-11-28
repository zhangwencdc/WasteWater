#!/usr/bin/perl
use strict;
#use warnings;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#Usage:基于bowite合并结果，提取耐药信息
####解读耐药报告
####中文注释部分

my %name;
#####仅包含了部分常见的抗生素英文翻中文
$name{"Amoxicillin"}="阿莫西林";
$name{"Piperacillin"}="哌拉西林";
$name{"Ampicillin"}="氨苄西林";
$name{"Ticarcillin"}="替卡西林";
$name{"Ceftazidime"}="头孢他啶";
$name{"Cefepime"}="头孢吡肟";
$name{"Cefotaxime"}="头孢噻肟";
$name{"Aztreonam"}="氨曲南";
$name{"Ceftriaxone"}="头孢曲松";
$name{"Ampicillin+Clavulanic acid"}="氨苄西林+克拉维酸";
$name{"Cefoxitin"}="头孢西丁";
$name{"Ertapenem"}="厄他培南";
$name{"Amoxicillin+Clavulanic acid"}="阿莫西林+克拉维酸";
$name{"Ticarcillin+Clavulanic acid"}="替卡西林+克拉维酸";
$name{"Sulfamethoxazole"}="复方新诺明";
$name{"Meropenem"}="美罗培南";
$name{"Imipenem"}="亚胺培南";
$name{"Piperacillin+Tazobactam"}="哌拉西林+他唑巴坦";
$name{"Trimethoprim"}="甲氧苄啶";
$name{"Ciprofloxacin"}="环丙沙星";
$name{"Tetracycline"}="四环素";
$name{"Doxycycline"}="多西环素";
$name{"Streptomycin"}="链霉素";
$name{"Fosfomycin"}="磷霉素";
$name{"Cephalotin"}="头孢噻酚Cephalotin";
$name{"Tobramycin"}="妥布霉素";
$name{"Gentamicin"}="庆大霉素";
$name{"Cephalothin"}="头孢噻酚Cephalothin";
$name{"Nalidixic acid"}="萘啶酸";
$name{"Colistin"}="粘菌素";
$name{"Tigecycline"}="替加环素";
$name{"Azithromycin"}="阿奇霉素";
$name{"Amikacin"}="阿米卡星";

#####
my $data=$Bin."/phenotypes.txt";

my $file=$ARGV[0];### Res.cat.90
my $sample=$ARGV[1];#Sample.info

#读sample信息
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

#读耐药信息
open(FILE,$data);
my %anno; my %tmp;my %class;my %class_n;
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	$class{$a[0]}=$a[1];
	$class_n{$a[1]}++;
	$anno{$a[0]}=$a[2];
	my $b=$anno{$a[0]};
	my @b=split",",$b;
	foreach my $bb (@b) {
		unless(substr($bb,0,1)=~/[0-9a-zA-Z]/){$bb=substr($bb,1);}
		unless(substr($bb,length($bb)-1,1)=~/[0-9a-zA-Z]/){$bb=substr($bb,0,length($bb)-1);}
		$tmp{$bb}++;
	}
	
}
close FILE;
open(F,$file);
my $out=$file.".detail";
open(OUT,">$out");my %drug;my %gene;my %gene_country;
while(1){
	my $line=<F>;
	unless($line){last;}
	chomp $line;
	my @a=split",",$line;
	my $n=basename($a[0]);
	my @b=split"\\.",$n;
	my $name=$b[0];
	unless(exists $genome{$name}){next;}
	my $class=$class{$a[1]};
	my $anno=$anno{$a[1]};
	print OUT  "$name,$a[1],$a[2],$a[3],$a[4],$a[5],Class: $class,Phenotype: $anno\n";
	$drug{$name}{$class}++;#print "$name,$class,$drug{$name}{$class}\n";
	$gene{$a[1]}++;
	my $country=$country{$name};$gene_country{$a[1]}{$country}++;
}
close F;
close OUT;

$out=$file.".stat";
open(OUT,">$out");

my @sample=sort keys %genome;
my @class=sort keys %class_n;my %sum;my %s_country;
foreach my $sample (@sample) {
	foreach my $class (@class) {
		#print "$sample,$class,$drug{$sample}{$class}\n";
		if(exists $drug{$sample}{$class}){
			print "$sample,$class,$drug{$sample}{$class}\n";
			$sum{$class}++;
			my $country=$country{$sample};
			$s_country{$class}{$country}++;
		}
	}
}
print OUT "Class\tSum\t";
my @country=sort keys %country_n;
foreach my $country (@country) {
	print OUT "$country\t";
}
print OUT "\n";


foreach my $class (@class) {
	print OUT "$class\t$sum{$class}\t";
	foreach my $country (@country) {
		print OUT "$s_country{$class}{$country}\t";
	}
	print OUT "\n";
}

my @gene=sort keys %class;
print OUT "Gene\tSum\t";

foreach my $country (@country) {
	print OUT "$country,";
}
print OUT "\n";


foreach my $gene (@gene) {
	print OUT "$gene,$gene{$gene},";
	foreach my $country (@country) {
		print OUT "$gene_country{$gene}{$country},";
	}
	print OUT "\n";
}