#!/usr/bin/perl
use strict;
#use warnings;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#Usage:����bowite�ϲ��������ȡ��ҩ��Ϣ
####�����ҩ����
####����ע�Ͳ���

my %name;
#####�������˲��ֳ����Ŀ�����Ӣ�ķ�����
$name{"Amoxicillin"}="��Ī����";
$name{"Piperacillin"}="��������";
$name{"Ampicillin"}="��������";
$name{"Ticarcillin"}="�濨����";
$name{"Ceftazidime"}="ͷ�����";
$name{"Cefepime"}="ͷ�����";
$name{"Cefotaxime"}="ͷ�����";
$name{"Aztreonam"}="������";
$name{"Ceftriaxone"}="ͷ������";
$name{"Ampicillin+Clavulanic acid"}="��������+����ά��";
$name{"Cefoxitin"}="ͷ������";
$name{"Ertapenem"}="��������";
$name{"Amoxicillin+Clavulanic acid"}="��Ī����+����ά��";
$name{"Ticarcillin+Clavulanic acid"}="�濨����+����ά��";
$name{"Sulfamethoxazole"}="������ŵ��";
$name{"Meropenem"}="��������";
$name{"Imipenem"}="�ǰ�����";
$name{"Piperacillin+Tazobactam"}="��������+�����̹";
$name{"Trimethoprim"}="�������";
$name{"Ciprofloxacin"}="����ɳ��";
$name{"Tetracycline"}="�Ļ���";
$name{"Doxycycline"}="��������";
$name{"Streptomycin"}="��ù��";
$name{"Fosfomycin"}="��ù��";
$name{"Cephalotin"}="ͷ�����Cephalotin";
$name{"Tobramycin"}="�ײ�ù��";
$name{"Gentamicin"}="���ù��";
$name{"Cephalothin"}="ͷ�����Cephalothin";
$name{"Nalidixic acid"}="�����";
$name{"Colistin"}="ճ����";
$name{"Tigecycline"}="��ӻ���";
$name{"Azithromycin"}="����ù��";
$name{"Amikacin"}="���׿���";

#####
my $data=$Bin."/phenotypes.txt";

my $file=$ARGV[0];### Res.cat.90
my $sample=$ARGV[1];#Sample.info

#��sample��Ϣ
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

#����ҩ��Ϣ
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