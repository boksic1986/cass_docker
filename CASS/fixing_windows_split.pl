#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw/mean median/;

use PerlIO::gzip;
use Getopt::Long;

my $usage="
	Options:
		-w <INT>		The windows size;
		-s <INT>		The windows slide size;
		-o <STR>		The output file;
		-ref <STR>		The reference fasta [/home/yxtkttkl/database/hg19/hg19_clean.fa]
		-h|?			Help!
";
my ($window_size,$slide_size,$out,$help);
my $ref="/home/yxtkttkl/database/hg19/hg19_clean.fa";
GetOptions(
	'w=i' => \$window_size,
	's=i' => \$slide_size,
	'o=s' => \$out,
	'ref=s' => \$ref,
	'h|?' => \$help
);

if($help or !$out){die "$usage\n";}
open REF,"<$ref";
open OUT,">:gzip","$out";
my ($chr,$seq);
while(<REF>){
	chomp;
	my $line=$_;
	if($line=~/^>/){
		$line=~s/>chr//ig;
		$line=$line=~/X/?23:$line=~/Y/?24:$line=~/M/?25:$line;
		if($.==1){$chr=$line;next;}
		my $chr_length=length $seq;
		for(my $i=1;$i<=$chr_length;$i+=$slide_size){
			next if $i+$window_size-1>$chr_length;
			my $window_seq=substr($seq,$i,$window_size);
			next if $window_seq=~/N/ig;
			my $window_gc=($window_seq=~s/[GC]//ig || 0);
			print OUT "$chr\t$i\t",$i+$window_size-1,"\t",$window_gc/$window_size,"\n";
		}
		undef $seq;
		$chr=$line;
	}
	$seq.=$line;
}
my $chr_length=length $seq;
for(my $i=1;$i<=$chr_length;$i+=$slide_size){
	next if $i+$window_size-1>$chr_length;
	my $window_seq=substr($seq,$i,$window_size);
	next if $window_seq=~/N/ig;
	my $window_gc=($window_seq=~s/[GC]//ig || 0);
	print OUT "$chr\t$i\t",$i+$window_size-1,"\t",$window_gc/$window_size,"\n";
}
close OUT;
close REF;
