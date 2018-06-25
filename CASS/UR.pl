#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw/mean median/;

use PerlIO::gzip;
use Getopt::Long;

my $usage="
	use to: calculate the ur for each sample;
	Options:
		-i/il <STR>				The bam file or list;
		-o <STR>				The output file;
		-t <STR>				The target region [whole chr];
		-h|?					Help!
";
my ($in,$in_list,$out,$target,$help);
GetOptions(
	'i=s' => \$in,
	'il=s' => \$in_list,
	'o=s' => \$out,
	't=s' => \$target,
	'h|?' => \$help
);

if($help or !$out){die "$usage\n";}
my %target;
if($target){
	open TAR,"<$target";
	while(<TAR>){
		chomp;
		my @line=split;
		$line[0]=~s/chr//g;
		$line[0]=($line[0] eq "X")?23:($line[0] eq "Y")?24:($line[0] eq "M")?25:$line[0];
		$target{$line[0]}{$line[1]}=$line[2];
	}
	close TAR;
}
if($in_list){
	system "mkdir -p $out" unless -d $out;
	open IL,"<$in_list";
	while(<IL>){
		chomp;
		my ($name,$in_path)=(split)[0,-1];
		unless(-f $in_path){next;}
		UR($in_path,"$out/$name.ur");
	}
	close IL;
}
if($in){UR($in,$out);}
sub UR{
	my ($ur_in,$ur_out)=@_;
	open IN,"samtools view $ur_in |";
	open OUT,">$ur_out";
	my ($total)=(0);
	my %ur;
	while(<IN>){
		chomp;
		my @line=split;
		my $line=$_;
		next if /^\@/;
		next if $line[5] eq "*";
		$line[2]=~s/chr//g;
		$line[2]=($line[2] eq "X")?23:($line[2] eq "Y")?24:($line[2] eq "M")?25:$line[2];
		next unless $line[2]=~/^\d+$/;
		next if $line[2]==25;
		if($line=~/XT\:A\:/){
			next unless $line=~/XT\:A\:U/;
		}else{next if $line=~/XA\:Z\:/;}
		$total++ if $line[2]<23;
		if($target){
			if(my ($start)=grep {$line[3]<=$target{$line[2]}{$_} and $line[3]>=$_} (keys %{$target{$line[2]}})){
				$ur{$line[2]}{$start}{$target{$line[2]}{$start}}++;
			}
		}else{$ur{$line[2]}++;}
	}
	foreach my $chr(sort {$a <=> $b} keys %ur){
		if($target){
			foreach my $s(sort {$a <=> $b} keys %{$ur{$chr}}){
				foreach my $e(keys %{$ur{$chr}{$s}}){
					print OUT "$chr\t$s\t$e\t$ur{$chr}{$s}{$e}\t$total\t",100*$ur{$chr}{$s}{$e}/$total,"\n";
				}
			}
		}else{print OUT "$chr\t$ur{$chr}\t$total\t",100*$ur{$chr}/$total,"\n";}
	}
	close IN;
	close OUT;
}
sub mean_sd{
	my @data=@_;
	my ($sum,$sd_sum)=(0,0);
	foreach my $d(@data){
		$sum+=$d;
	}
	my $aver=$sum/@data;
	foreach my $dd(@data){
		$sd_sum+=($dd-$aver)**2;
	}
	my $sd=sqrt($sd_sum/@data-1);
	return $aver,$sd;
}
