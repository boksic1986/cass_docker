#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;

my $usage="
	Usage:	perl ext.pl -i <IN sorted and rmdupped bam file> -o <*.ext>
	Motified Date: 2016-07-19
	Options:
		-i <STR>		The format of input file is *.bam(BWA/Samtools);
		-o <STR>		The output file (*.ext)
		-gender_cutoff <STR>	The cutoff of gender identification [0.125];
		-h|?			Help!
";
my ($in,$out,$help);
my $gender=0.125;
GetOptions(
	'i=s' => \$in,
	'o=s' => \$out,
	'gender_cutoff=s' => \$gender,
	'h|?' => \$help
);
if($help or !$in){die "$usage\n";}
extract();
sub extract{
	open IN,"samtools view $in |";
	open OUT,">$out";
#	my ($total,$chrY)=(0,0);
	while(<IN>){
		chomp;
		my @line=split;
		my $line=$_;
		next if /^\@/;
		next if $line[2] eq "*";
		if($line=~/\s+XT\:A\:/){
			if($line!~/\s+XT\:A\:U/){next;}
		}elsif($line=~/\s+XA\:Z\:/){next;}
		$line[2]=~s/chr//g;
		$line[2]=($line[2] eq "X")?23:($line[2] eq "Y")?24:($line[2] eq "M")?25:$line[2];
		next unless $line[2]=~/^\d+$/;
		next if $line[2]==25;
#		$total++;
#		$chrY++ if $line[2]==24;
		my $length=length $line[9];
		my $gc=($line[9]=~s/[GC]//ig || 0);
		print OUT "$line[2]\t$line[3]\t$gc\t$length\n";
	}
#	print $chrY/$total>$gender?"M":"F";
	close IN;close OUT;
	system "sort -n -k 1 -k 2 -S 8g $out > $out.sort";
	system "mv $out.sort $out";
	system "gzip $out";
}
