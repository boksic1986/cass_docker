#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use Cwd qw/abs_path/;

my $usage="
	Options:
		-contig <STR>		The contig file;
		-i <STR>		The input file list (fastq files);
		-o <STR>		The output folder;
		-bwa			The input files are alignment results [fastq];
		-h|?			Help!
";
my ($contig,$in,$out,$bwa_flag,$help);
GetOptions(
	'i=s' => \$in,
	'o=s' => \$out,
	'contig=s' => \$contig,
	'bwa' => \$bwa_flag,
	'h|?' => \$help
);

if($help or !$contig){die "$usage\n";}
$out=$out?abs_path($out):abs_path();
system "mkdir -p $out" unless -d $out;
open IN,"<$in";
open CON,"<$contig";
open SHL,">$out/shell.sh";
#open CNVL,">$out/cnv.list";
open BAML,">$out/bam.list" unless $bwa_flag;
#open DATA,">$out/data.list";
## Read the configuration file;
my (%config);
while(<CON>){
	chomp;
	next if /^#/;
	my @line=split /=/,$_;
	$config{$line[0]}=$line[1] if $line[1];
}
## Read the input files;
while(<IN>){
	chomp;
	my @line=split;
	my $outputdir="$out/$line[0]";
	system "mkdir -p $outputdir" unless -d $outputdir;
	system "mkdir -p $outputdir/graph";
	system "mkdir -p $outputdir/graph/copyratio";
#	system "mkdir -p $outputdir/graph/copynumber";
	my $shfile="$outputdir/$line[0].sh";
	print SHL "nohup sh $shfile > $shfile.olog 2> $shfile.elog &\n";
	open SH,">$shfile";
	my $bwa=(split /\s+/,$config{"bwa"})[0];
#	print SH "perlbrew switch perl-5.14.4\n";
	if(!$bwa_flag){
		print SH $config{"cut"}," $line[1] $outputdir/$line[0].cut.fastq.gz\n" if (exists $config{"cut"});
		print SH $config{"cut"}," $line[-1] $outputdir/$line[0].R2.cut.fastq.gz\n" if (exists $config{"cut"} and @line>2);
		my $fastq1=(exists $config{"cut"})?"$outputdir/$line[0].cut.fastq.gz":$line[1];
		my $fastq2=(exists $config{"cut"})?"$outputdir/$line[0].R2.cut.fastq.gz":$line[-1] if @line>2;
		print SH $config{"bwa"}," $fastq1 -f $outputdir/$line[0].sai\n";
		print SH $config{"bwa"}," $fastq2 -f $outputdir/$line[0].R2.sai\n" if @line>2;
		print SH "$bwa ",@line>2?"sampe ":"samse ","-r \"\@RG\\tID:$line[0]\\tSM:$line[0]\\tPL:illumina\" ",$config{"reference"}," $outputdir/$line[0].sai ",@line>2?"$outputdir/$line[0].R2.sai ":"",$fastq1,@line>2?" $fastq2":""," -f $outputdir/$line[0].sam\n";
		print SH $config{"samtools"}," view -Sb $outputdir/$line[0].sam | ",$config{"samtools"}," sort -m 8g -T $outputdir > $outputdir/$line[0].bam\n";
		print SH $config{"samtools"}," rmdup $outputdir/$line[0].bam $outputdir/$line[0].rmdup.bam 2> $outputdir/$line[0].rmdup.stat\n";
		print SH $config{"extract"}," -i $outputdir/$line[0].rmdup.bam -o $outputdir/$line[0].ext\n";
	}else{
		print SH $config{"samtools"}," rmdup $line[-1] $outputdir/$line[0].rmdup.bam 2> $outputdir/$line[0].rmdup.stat\n";
		print SH $config{"extract"}," -i $outputdir/$line[0].rmdup.bam -o $outputdir/$line[0].ext\n";
	}
	print SH "gender=`",$config{"tags"}," -w ",$config{"windows"}," -i $outputdir/$line[0].ext.gz -o $outputdir/$line[0].tagsinfo.gz`\n";
	print SH $config{"seg"}," -n ",$config{"n_region"}," -band ",$config{"cytoBand"}," -i $outputdir/$line[0].tagsinfo.gz -o $outputdir/$line[0].copyratio.cnv.txt -g \$gender\n";
#	print SH $config{"seg"}," -n ",$config{"n_region"}," -band ",$config{"cytoBand"}," -i $outputdir/$line[0].tagsinfo.gz.CN.gz -o $outputdir/$line[0].copynumber.cnv.txt -fa 1.5 -fb 2.5 -g F\n";
	print SH $config{"drawcnv"}," -r ",$config{"R"}," -i $outputdir/$line[0].copyratio.cnv.txt -o $outputdir/graph/copyratio -z $outputdir/$line[0].tagsinfo.gz -ylim 0,3 -name $line[0]\n";
#	print SH $config{"drawcnv"}," -r ",$config{"R"}," -i $outputdir/$line[0].copynumber.cnv.txt -o $outputdir/graph/copynumber -z $outputdir/$line[0].tagsinfo.gz.CN.gz -ylim 0,6 -name $line[0]\n";
	print SH $config{"boxplot"}," -name $line[0] -i $outputdir/$line[0].tagsinfo.gz -o $outputdir/$line[0].chrratio.xls -g $outputdir/$line[0].chrratio.pdf -R ",$config{"R"},"\n";
	print SH $config{"karyomap"}," -i $outputdir/$line[0].copyratio.cnv.txt -o $outputdir/$line[0].cnv.svg -gender \$gender\nconvert $outputdir/$line[0].cnv.svg $outputdir/$line[0].cnv.png\n";
#	print SH $config{"samtools"}," view $outputdir/$line[0].bam |awk '{print \$3,\$10}'|grep \"chr\"|awk '{print \$2}'|sort|uniq -u|wc -l >> $outputdir/$line[0].data\n";
	print SH "rm $outputdir/$line[0].sam $outputdir/$line[0].sai $outputdir/$line[0].ext.gz\n";
	close SH;
	#print CNVL "$line[0]\t$outputdir/$line[0].cnv\n";
	print BAML "$line[0]\t$outputdir/$line[0].bam\n" unless $bwa_flag;
#	print DATA "$line[0]\t$outputdir/$line[0].data\n";
}
## QC & make data production;
#open QC,">$out/QCsummary.sh";
#print QC $config{"QCsummary"}," -il $out/data.list -o $out/dataproduction.xls\n";
## Run the shells;
close IN;close CON;close SHL;
#close CNVL;
close BAML unless $bwa_flag;
#close QC;
