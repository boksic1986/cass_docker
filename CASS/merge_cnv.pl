#!/usr/bin/perl -w

## For: CNV filtering;
## Author: Chunlei Zhang;
## Email: yxtkttkl@gmail.com;
## Version: v17.8

use strict;
use Statistics::Distributions;
use Statistics::Sequences::Runs;
use Statistics::Data;
use Getopt::Long;
use List::MoreUtils qw/uniq/;
use List::Util qw/sum min max/;

my $usage=<<'USAGE';
Usage: 

	-cnv <STR>		The raw CNV file;
	-tags <STR>		The tags file [false];
	-size <INT>		The cutoff of CNV size [5000000];
	-o <STR>		The output file;
	-h | help | ?		Help!

USAGE

my ($cnv_input,$tags,$out,$help);
my $cnv_size=5000000;
GetOptions(
	'cnv=s' => \$cnv_input,
	'tags=s' => \$tags,
	'size=i' => \$cnv_size,
	'o=s' => \$out,
	'h|help|?' => \$help
);

die "$usage\n" if ($help || !$cnv_input || !$out);
open CNV,"<$cnv_input";
if($tags){if(-B $tags){open TAGS,"<:gzip","$tags";}else{open TAGS,"<$tags";}}
open OUT,">$out";
chomp(my $title=<CNV>);
my @cnv;
while(<CNV>){
	chomp;
	my $line=$_;
	push @cnv,$line;
}
for(my $i=1;$i<$#cnv;$i++){
	my @last_line=split /\s+/,$cnv[$i-1];
	my @line=split /\s+/,$cnv[$i];
	if($line[0]==$last_line[0]){
		if($line[5] eq $last_line[5]){
			$line[1]=$last_line[1];
			
		}
	}else{
		print OUT $last_line,"\n";
		next;
	}
}
