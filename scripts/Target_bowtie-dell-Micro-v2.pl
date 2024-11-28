#!/use/bin/perl

=head1 Name: virus鉴定流程

=head1 Description
可用于指定基因or基因组序列的识别
输出文件：鉴别为该基因的序列数量及fasta格式的序列。
=head1

  Author: zhangwen, zhangwen@icdc.cn
  Version: 1.0, Date: 2019-5-28

=head1 Usage:

perl  Target_bowtie.pl -1 fq1.fastq -2 fq2.fastq -o outputfile
perl  Target_bowtie.pl -1 fq1.fastq -2 fq2.fastq -o outputfile -m 100 -c 80
=head1 Example
=cut

use strict;
#use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
#use warnings;
use Pod::Text;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;
##get options from command line into variables and set default values
###-----------------Time record-----------------------------###
my $Time_Start = sub_format_datetime(localtime(time())); #......
my $Data_Vision = substr($Time_Start,0,10);

###针对不可拼接的双向测序数据
my $fq1;my $fq2;my $o;my $match_cutoff;my $cutoff;my $target;
GetOptions("1|fq1:s" => \$fq1,
"2|fq2:s" => \$fq2,
"T|Target:s" => \$target,
"m|ML:s" => \$match_cutoff,
"c|C:s" => \$cutoff,
"o|out:s" =>\$o
);
die `pod2text $0` if ( !defined $fq1 || !defined $target);
if ( !defined $o ){$o="./test";};
###与virus的比对参数 可调整
if ( !defined $match_cutoff ){ $match_cutoff=70;}###最小比对长度70
if ( !defined $cutoff ){ $cutoff=80; }##最小比对百分百80%


###bowtie比对库路径  程序迁移后需相应修改
my $blast="blastn";
my $bowtie2build="/home/zhangwen/bin/bowtie2-2.4.4-linux-x86_64/bowtie2-build";
my $bowite="/home/zhangwen/bin/bowtie2-2.4.4-linux-x86_64/bowtie2";

unless(-e "$target.1.bt2"){system "$bowtie2build $target $target\n";}

open(FILE,"$target");
my %gene;my %genelen;my $len;my $name;
while(1){
	my $line=<FILE>;
	unless($line){$genelen{$name}=$len;#print "$name,$len\n";
	last;}
	chomp $line;
	unless(substr($line,length($line)-1,1)=~/[0-9a-zA-Z]/){$line=substr($line,0,length($line)-1);}
	unless($line=~/>/){$len+=length($line);next;}
	$genelen{$name}=$len;
	#print "$name,$len\n";
	my @a=split" ",$line;
	$name=substr($a[0],1);$gene{$name}=$line;
	
	
	$len=0;
}
close FILE;


###virus 比对
print "Step1: Target align\n";
if(defined $fq2){system "$bowite -1 $fq1 -2 $fq2 -x $target -S $o.target.sam --no-unal\n";}else{system "$bowite -U $fq1  -x $target -S $o.target.sam --no-unal\n";}
system "samtools view -bS $o.target.sam > $o.target.bam\n";
system "samtools sort $o.target.bam -o $o.target.sort.bam\n";
system "samtools index $o.target.sort.bam\n";
system "bcftools mpileup -Ou -f  $target $o.target.sort.bam | bcftools call -mv -Ov -o  $o.target.vcf\n";

###
open(FILE,"$o.target.sam");
#print "$fq1,$fq2,$o\n";
#print "$o.virus.sam\n";
my %read;
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	if(substr($line,0,1) eq "@"){next;}
	my @a=split" ",$line;
	unless($a[2]=~/[0-9a-zA-Z]/){next;}
	my $len=length($a[9]);
	
	my $tmp;
	foreach my $a (@a) {
		if($a=~/^MD:Z:/){$tmp=$a;}
	}
	my $l=length($tmp);
	my $num="";my $match=0;my $total=$len;
	foreach  (0..($l-1)) {
		my $site=substr($tmp,$_,1);
		if($site=~/[0-9]/){
			$num=$num.$site;
		}else{
			$match+=$num;

			$num="";			
		}
	}
	$match+=$num;
	
	unless($len>0){next;}
	unless($match>0){next;}
	my $maper=$match/$len*100;
	if($match<$match_cutoff && $maper<$cutoff){next;}  ###
	my @b=split"/",$a[0];
	$read{$b[0]}++;


}
close FILE;


my @key=keys %read; my %virus;
my %result; my $num=0;
open(OL,">$o.matchread");
foreach my $key (@key) {
	if($read{$key}<1){next;} #没有与virus比对
	
	$result{$key}++;
	print OL "$key\n";
	$num++;
}
close OL;
print "Candidate paired reads uniq map to Target is $num\n";

   open(OUT,">$o.fq1.filter.fasta"); #open(O1,">$o.fq1.filter.fq");
open(FILE,$fq1);my $abs;
while(1){
        my $line=<FILE>;
        unless($line){last;}
        chomp $line;
        my $seq=<FILE>;
        chomp $seq;
        my $a=<FILE>;
        chomp $a;
        my $b=<FILE>;
        chomp $b;
        $line=substr($line,1);
        my @name=split" ",$line;
        my @b=split"/",$name[0];
        my $name=$b[0];

      # print "$name,$read{$name}\n";
		$abs++;
        if(exists $read{$name}){
               print OUT ">";
                print OUT "$name fq1\n$seq\n";
#				print O1 "@";
#				print O1 "$name\n$seq\n$a\n$b\n";
        }

}

close FILE;close OUT;

if(defined $fq2){
	open(OUT,">$o.fq2.filter.fasta");#open(O2,">$o.fq2.filter.fq");
open(FILE,$fq2);
while(1){
        my $line=<FILE>;
        unless($line){last;}
        chomp $line;
        my $seq=<FILE>;
        chomp $seq;
        my $a=<FILE>;
        chomp $a;
        my $b=<FILE>;
        chomp $b;
        $line=substr($line,1);
        my @name=split" ",$line;
        my @b=split"/",$name[0];
        my $name=$b[0];

#       print "$name\n";
        if(exists $result{$name}){
			print OUT ">";
                print OUT "$name fq2\n$seq\n";
#				print O2 "@";
#				print O2 "$name\n$seq\n$a\n$b\n";
        }

}
close FILE;
}
close OUT;
print "Step2: Target Gene filter align\n";
if(defined $fq2){system "$bowite -1 $o.fq1.filter.fasta -2 $o.fq2.filter.fasta -x $target -S $o.target.filter.sam -f --no-unal\n";}else{system "$bowite -U $o.fq1.filter.fasta  -x $target -S $o.target.filter.sam -f --no-unal\n";}
open(FILE,"$o.target.filter.sam");
#print "$fq1,$fq2,$o\n";
#print "$o.virus.sam\n";
my %read;my %start;my %end;my %align;
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	if(substr($line,0,1) eq "@"){next;}
	my @a=split" ",$line;
	unless($a[2]=~/[0-9a-zA-Z]/){next;}
	my $len=length($a[9]);
	
	my $tmp;
	foreach my $a (@a) {
		if($a=~/^MD:Z:/){$tmp=$a;}
	}
	my $l=length($tmp);
	my $num="";my $match=0;my $total=$len;
	foreach  (0..($l-1)) {
		my $site=substr($tmp,$_,1);
		if($site=~/[0-9]/){
			$num=$num.$site;
		}else{
			$match+=$num;

			$num="";			
		}
	}
	$match+=$num;
	
	unless($len>0){next;}
	unless($match>0){next;}
	my $maper=$match/$len*100;
	if($match<$match_cutoff && $maper<$cutoff){next;}  ###
	$read{$a[2]}++;
	if($a[1] == 16){ ###与负链匹配
			foreach  (($a[3]-length($a[9])+1)..$a[3]) {
				$align{$a[2]}{$_}++;
		}
	
	}else{
		
		foreach  ($a[3]..($a[3]+length($a[9])-1)) {
				$align{$a[2]}{$_}++;
		}
	}
	

}
close FILE;


open(OUT,">$o.bowtie.report");
print OUT "Gene ID,Target Length,Match Read Num,Align Length,Coverage%,Depth,Gene Annotation,\n";
my @key=sort {$read{$b}<=>$read{$a}} keys %read;
foreach my $key (@key) {
	my $align;my $depth;
	foreach  (1..$genelen{$key}) {
		if($align{$key}{$_}>=1){$align++;$depth+=$align{$key}{$_};}
	}
	if($align>0){$depth=$depth/$align;}
	unless($genelen{$key}>0){next;}
	my $coverage=$align/$genelen{$key}*100;
	$coverage=sprintf "%.2f",$coverage;
	
	print OUT "$key,$genelen{$key},$read{$key},";
	print OUT "$align,$coverage %,$depth,$gene{$key}\n";
}
close OUT;

my $Time_End= sub_format_datetime(localtime(time()));
print "Running from [$Time_Start] to [$Time_End]\n";



#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+



sub sub_format_datetime #.....
{
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
        $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon, $day, $hour, $min, $sec);
}
#####
sub gc{
        my $seq=shift @_;
        my $gc=$seq=~tr/(G|C|g|c)/(G|C|g|c)/;
        my $l=length($seq);


        return ($gc,$l);
}

