#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;

my $usage="
	Usage:	perl Sta_ext.pl -i <IN sorted bam file>
	Output:	*.ext (GC info ) && STAT_*.txt (Total info)
	Options:
		-i				The format of input file is *.bam(BWA/Samtools);
		-h|?			Help!
";
my ($in,$help);
GetOptions(
	'i=s' => \$in,
	'h|?' => \$help
);
if($help or !$in){die "$usage\n";}
my $out = $in;
$out =~s/.bam//g;

open OUT,">$out.ext";
open STAT,">$out.data";
Stat($in);
extract($in);
my $dir=dirname($in);
`sort -n -k 1 -k 2 -S 1g -T $dir $out.ext > $out.sort`;
`mv $out.sort $out.ext`;
system "gzip $out.ext";
sub extract{
	my ($in)=@_;
	open IN,"/usr/bin/samtools view $in |";
	my $line="";
	my @line ="";
	my $record ="";
	my @record ="";	
	while(<IN>){
		chomp;
		my @line=split;
		my $line=$_;
		next if /^\@/;
		next if $line[5] eq "*";
		next unless $line=~/XT\:A\:U/;
		$line[2]=~s/chr//g;
		$line[2]=($line[2] eq "X")?23:($line[2] eq "Y")?24:($line[2] eq "M")?25:$line[2];
		next unless $line[2]=~/^\d+$/;
		next if $line[2]==25;		
		my $length=length $line[9];
		my $gc=($line[9]=~s/[GC]//ig || 0);
		$record =$line[2]."\t".$line[3]."\t".$gc."\t".$length;
		push(@record,$record);
	}
	my @unique = uniq( @record );
	sub uniq {
		my %seen;
		return grep { !$seen{$_}++ } @_;
	}
	my $CHR ="";
	my $POS =0;
	my $GC =0;
	my $LEN =0;	
	my $Avalible_reads_num =0;
	my $GC_Avalible = 0;
	my $Base_len_Available =0;
	my $GC_Avail_PER =0;
	foreach my $unique(@unique){
		next if $unique =~ /^\s*$/;
		print OUT "$unique\n";
		$Avalible_reads_num++;
		my ($CHR,$POS,$GC,$LEN)=split("\t",$unique);
		$GC_Avalible = $GC_Avalible+$GC;
		$Base_len_Available = $Base_len_Available+$LEN;	
	}
	close IN;
	$GC_Avail_PER = $GC_Avalible/$Base_len_Available;
	print STAT $Avalible_reads_num."\n";
	print STAT $GC_Avail_PER."\n";
}

sub Stat{
	my ($in)=@_;
	open IN,"/usr/bin/samtools view $in |";
	my (@bases,@line);
	my $line = "";
	my $Base_count=0;
	my $bases = "";
	my $temp =0;
	my $Q20=0;
	my $Q30=0;
	my $GC_temp = "";
	my $seq_len =0;
	my $gc=0;
	my $Total_Reads_num =0 ;
	my $Algn_reads_num =0;
	my $Q20_PER =0;
	my $Q30_PER =0;
	my $Base_len =0;
	my $Base_len_aln =0;
	my $GC_total=0;
	my $GC_total_aln =0;
	my $GC_aln_PER=0;
	my $GC_total_PER=0;
	while(<IN>){
		chomp;
		$Total_Reads_num ++;
		@line=split;
		$line=$_;
		@bases=split//,$line[10];
		foreach $bases(@bases){
			$Base_count ++;
			$temp = ord($bases)-33;	
			if ($temp >= 20){
				$Q20++ ;
			}
			if($temp >= 30)	{
				$Q30++;
			}
		}
		##Get GC%
		$GC_temp = $line[9];
		$seq_len = length($GC_temp);
		$gc= ($GC_temp=~s/[GC]//g || 0);
		$Base_len = $Base_len + $seq_len;
		$GC_total = $GC_total + $gc;
		if ($line[2] =~ /chr/){
			$Algn_reads_num ++;
			$GC_temp = $line[9];
			$seq_len = length($GC_temp);
			$gc= ($GC_temp=~s/[GC]//g || 0);
			$Base_len_aln = $Base_len_aln + $seq_len;
			$GC_total_aln = $GC_total_aln + $gc;	
	}
	}

	$Q20_PER = $Q20/$Base_count;
	$Q30_PER = $Q30/$Base_count;
	$GC_total_PER = $GC_total/$Base_len;
	$GC_aln_PER = $GC_total_aln/$Base_len_aln;
	print STAT "$Total_Reads_num\n";
	print STAT "$Algn_reads_num\n";
	print STAT "$Q20_PER\n";
	print STAT "$Q30_PER\n";
	print STAT "$GC_total_PER\n";
	print STAT "$GC_aln_PER\n";
	close IN;
}
		
