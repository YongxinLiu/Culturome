#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

# 编写分菌鉴定mappingfile(metadata)的脚本
# 通常每个文库有48个板index，每板96孔
# 基本信息为#SampleID       BarcodeSequence LinkerPrimerSequence    ReversePrimer   barcodeF        barcodeR        plate
# 可选信息有variety	compartment	medium    batch   species Description

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq!
Usage:    write_mapping_file.pl -i ForwardBarcodeFile -b ReverseBarcodeFile -L LibraryID -o Metadata -F ForwardPrimer -R RevsersePrimer -p PlateNumber -v Variety -c Compartment -m Medium -B batch -s Species
Function: write a mapping file for identifying culture microbe
Command:  -i Forward barcode file for wells, default "script/barcodeF96.txt"
          -b Reverse barcode file file for plates, default "script/barcodeR48.txt"
          -L library name, default L1, such as L1, L2
          -o output metadata/mapping file name (Must), such as L1.txt, L2.txt
          -F forward primer, default 16S rDNA 799F "AACMGGATTAGATACCCKG" 
          -R reverse primer, default 16S rDNA 1193R "ACGTCATCCCCACCTTCC" 
          -p plate number, default or max is 48, range 1-48
          -v variety for plant, default Nippobare
          -c compartment type, default Root, such as Rhizosphere, Soil, Leaf and so on
          -m medium type, default TSB, such as R2A
          -B batch type, default 1, such as 1,2,3
          -s species, default Rice
          -d description, default WildType
Author:   Liu Yong-Xin, E-mail: yxliu\@genetics.ac.cn, wechat: yongxinliu, QQ: 42789409
Version:  v2.1
Update:   2020/3/23
Notes:    Change hole 1-96 to A1-H12, Add variety and medium between compartment column
Example: 
# Using all default parameters and write metadata "L1.txt"
write_mapping_file.pl -i script/barcodeF96.txt -b script/barcodeR48.txt -L L1 -o output/L1.txt
# The full paramters
write_mapping_file.pl -i script/barcodeF96.txt -b script/barcodeR48.txt \
    -L L1 -o input/L1.txt \
    -F AACMGGATTAGATACCCKG -R ACGTCATCCCCACCTTCC \
    -p 48 -v Nippobare -c Root -m TSB -B 1 -s Rice -d WildType
# -h for help and more details
write_mapping_file.pl -h
\n!
    )
}

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:b:p:F:R:c:s:B:L:P:v:m:', \%opts );
&usage unless ( exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{i}="script/barcodeF96.txt" unless defined($opts{i});
$opts{b}="script/barcodeR48.txt" unless defined($opts{b});
$opts{L}="L1" unless defined($opts{L});
$opts{o}="L1.txt" unless defined($opts{o});
$opts{F}="AACMGGATTAGATACCCKG" unless defined($opts{F});
$opts{R}="ACGTCATCCCCACCTTCC" unless defined($opts{R});
$opts{p}=48 unless defined($opts{p});
$opts{v}="Nippobare" unless defined($opts{v});
$opts{c}="Root" unless defined($opts{c});
$opts{m}="TSB" unless defined($opts{m});
$opts{B}=1 unless defined($opts{B});
$opts{s}="Rice" unless defined($opts{s});
$opts{d}="WildType" unless defined($opts{d});



# 读取正向Barcode列表
open INPUT,"<$opts{i}";
my @barcodeF; #database in array
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	push @barcodeF,$tmp[0];
}
close INPUT;

# 读取反向Barcode列表
open DATABASE,"<$opts{b}";
my @barcodeR; #database in array
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @barcodeR,$tmp[0];
}
close DATABASE;

# 生成96孔板编号
my @hole;
@line=("A","B","C","D","E","F","G","H");
my $i=1;
foreach  $l(@line) {
	foreach  $r(1..12) {
          $hole[$i]="$l$r";
#          print $hole[$i],"\n";
          $i++;
	}
}

###############################################################################
#Main text.
###############################################################################
# 输出4608行的mapping file
open OUTPUT,">$opts{o}";
print OUTPUT "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tbarcodeF\tbarcodeR\tlibrary\tplate\twell\tvariety\tcompartment\tmedium\tbatch\tspecies\tDescription\n";
foreach $i(1..$opts{p}) {
	foreach $j(1..96) {
          $F=$j-1;
          $R=$i-1;
          printf OUTPUT "$opts{L}P%02d$hole[$j]\t$barcodeF[$F]$barcodeR[$R]\t$opts{F}\t$opts{R}\t$barcodeF[$F]\t$barcodeR[$R]\t$opts{L}\tP%02d\t$hole[$j]\t$opts{v}\t$opts{c}\t$opts{m}\t$opts{B}\t$opts{s}\t$opts{d}\n",$i,$i;
	}
}
close OUTPUT;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

