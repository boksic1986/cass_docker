#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Cwd qw/abs_path/;
use Statistics::Basic qw (:all);
use PerlIO::gzip;
use Getopt::Long;

my $usage="
	Options:
		-i/il <STR>		The sorted input file(s) [bam/sam];
		-o <STR>		The output file;
		-g <STR>		The cutoff of gender identification [0];
		-h|?			Help!
	Example: perl $0 -il <bam.list> -o <dataproduction>
";
my ($in,$in_list,$out,$help);
my $samtools="/usr/bin/samtools";
my $gender_cutoff=0;
GetOptions(
	'i=s' => \$in,
	'il=s' => \$in_list,
	'o=s' => \$out,
	'g=s' => \$gender_cutoff,
	'h|?' => \$help
);

if($help or !$out){die "$usage\n";}
if($in){
	depth($in,$out);
}
if($in_list){
	open IL,"<$in_list";
	while(<IL>){
		chomp;
		my @line=split;
		my $file_out="$out/$line[0].DepthCov.txt";
		depth($line[-1],$file_out);
	}
	close IL;
}
sub depth{
	my ($din,$dout)=@_;
	my @chr_length=(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566,16571);
	if(-B $din){open IN,"samtools view $din |";}else{open IN,"<$din";}
	open OUT,">$dout";
	my (%chr_depth,%chr_cov,%chr_reads,$Y_cov);
	my ($total_depth,$total_cov,$total_chr_length,$last_end,$last_chr,$total_reads)=(0,0,0,0,"chr0",0);
	while(<IN>){
		chomp;
		my @line=split;
		my $line=$_;
		next if /^\@/;
		next if $line[2] eq "*";
		$line[2]=~s/chr//ig;
		$line[2]=$line[2]=~/X/?23:$line[2]=~/Y/?24:$line[2]=~/M/?25:$line[2];
		if($line=~/\s+XT\:A\:/){
			if($line!~/\s+XT\:A\:U/){next;}
		}elsif($line=~/\s+XA\:Z\:/){next;}
		$total_reads++ unless $line[2]==25;
		$chr_reads{$line[2]}++;
		my @ignore_locus=$line[5]=~/\d+S/g;
		@ignore_locus=map {/(\d+)/} @ignore_locus;
		my $end=@ignore_locus==2?($line[3]+length($line[9])-$ignore_locus[0]-$ignore_locus[1]-1):@ignore_locus==1?($line[3]+length($line[9])-$ignore_locus[0]-1):($line[3]+length($line[9])-1);
		$chr_depth{$line[2]}+=$end-$line[3]+1;
		$chr_cov{$line[2]}+=($line[3]>$last_end)?($end-$line[3]+1):($end>$last_end)?($end-$last_end):0;
		$last_end=($line[2] ne $last_chr)?0:$end;
		$last_chr=$line[2];
	}
	for(my $chr=1;$chr<=25;$chr++){
		my $chr_length=$chr_length[$chr-1];
		my $aver_depth=(exists $chr_depth{$chr})?($chr_depth{$chr}/$chr_length):0;
		$total_depth+=(exists $chr_depth{$chr})?$chr_depth{$chr}:0;
		my $aver_cov=(exists $chr_cov{$chr})?($chr_cov{$chr}/$chr_length):0;
		$total_cov+=(exists $chr_cov{$chr})?$chr_cov{$chr}:0;
		my $aver_depth_coveraged=(exists $chr_depth{$chr})?($chr_depth{$chr}/$chr_cov{$chr}):0;
		$Y_cov=$aver_cov if $chr==24;
		$total_chr_length+=$chr_length;
		my $ur=(exists $chr_reads{$chr})?(100*$chr_reads{$chr}/$total_reads):0;
		my $chr_name=$chr;
		$chr_name=$chr_name==23?"X":$chr_name==24?"Y":$chr_name==25?"M":$chr_name;
		print OUT "chr$chr_name:$chr_length\tAverage_Depth:$aver_depth\tCoveraged_Depth:$aver_depth_coveraged\tCoverage:$aver_cov\tUR%:$ur\n";
	}
	print OUT "genome:$total_chr_length\tAverage_Depth:",$total_depth/$total_chr_length,"\tCoveraged_Depth:",$total_depth/$total_cov,"\tCoverage:",$total_cov/$total_chr_length,"\n";
	print OUT "Gender:\t",$Y_cov>$gender_cutoff?"M":"F","\n";
	close IN;
	close OUT;
}
