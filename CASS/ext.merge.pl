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
		-i <STR>		The input file list, the format is *.bam(BWA/Samtools);
		-o <STR>		The output file (*.ext)
		-samtools <STR>		The samtools path [/usr/bin/samtools];
		-h|?			Help!
";
my ($in,$out,$help);
my $samtools="/usr/bin/samtools";
GetOptions(
	'i=s' => \$in,
	'o=s' => \$out,
	'samtools=s' => \$samtools,
	'h|?' => \$help
);
if($help or !$in){die "$usage\n";}
extract();
sub extract{
	open IL,"<$in";
	open OUT,">$out";
	my $ext_line_num=0;
	while(<IL>){
		chomp;
		my @line=split;
		open IN,"$samtools view $line[-1] |";
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
			my $length=length $line[9];
			my $gc=($line[9]=~s/[GC]//ig || 0);
			print OUT "$line[2]\t$line[3]\t$gc\t$length\n";
			$ext_line_num++;
		}
		close IN;
	}
	close OUT;
	close IL;
	print $ext_line_num,"\n";
	system "sort -n -k 1 -k 2 -S 8g $out > $out.sort";
#	system "mv $out.sort $out";
#	system "gzip $out";
}
