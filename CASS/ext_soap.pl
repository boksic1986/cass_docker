#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;

my $usage="
	Usage:	perl ext.pl -i <IN sorted and rmdupped soap file> -o <*.ext>
	Motified Date: 2016-07-19
	Options:
		-i <STR>		The format of input file is *soap;
		-o <STR>		The output file (*.ext)
		-h|?			Help!
";
my ($in,$out,$help);
GetOptions(
	'i=s' => \$in,
	'o=s' => \$out,
	'h|?' => \$help
);
if($help or !$in){die "$usage\n";}
extract();
sub extract{
	if(-B $in){open IN,"<:gzip","$in";}else{open IN,"<$in";}
	open OUT,">$out";
	while(<IN>){
		chomp;
		my @line=split;
		next unless $line[3]==1;
		$line[7]=~s/chr//g;
		$line[7]=($line[2] eq "X")?23:($line[2] eq "Y")?24:($line[2] eq "M")?25:$line[2];
		next unless $line[7]=~/^\d+$/;
		next if $line[7]==25;
		my $length=length $line[1];
		my $gc=($line[1]=~s/[GC]//ig || 0);
		print OUT "$line[7]\t$line[8]\t$gc\t$length\n";
	}
	close IN;close OUT;
	system "sort -n -k 1 -k 2 -S 8g $out > $out.sort";
	system "mv $out.sort $out";
	system "gzip $out";
}
