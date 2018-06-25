#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 <cut size> <fastq.gz> <outputfile.gz>\n" unless @ARGV;
if(-B $ARGV[1]){open IN,"<:gzip","$ARGV[1]";}else{open IN,"<$ARGV[1]";}
open OUT,">:gzip","$ARGV[2]";
while(<IN>){
	chomp;
	my $line=$_;
	if($.%2==0){
		my $new_line=substr($line,$ARGV[0]);
		print OUT "$new_line\n";
	}else{
		print OUT "$line\n";
	}
}
close IN;close OUT;
