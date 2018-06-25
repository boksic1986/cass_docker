#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;
use Getopt::Long;

my $usage="
	Options:
		-il <STR>				The data.list;
		-o <STR>				The output file;
		-h|?					Help!
";
my ($in_list,$out,$help);
GetOptions(
        'il=s' => \$in_list,
        'o=s' => \$out,
        'h|?' => \$help
);

if($help or !$in_list){die "$usage\n";}
open IL,"<$in_list";
open OUT,">$out";
print OUT "Sample\tQ20\tQ30\tGC(fastq)\tGC(aln)\t#Reads\tAlignment%\tduplication\n";
while(<IL>){
	chomp;
	my @line=split;
	open IN,"<$line[-1]";
	chomp(my @data=<IN>);
	close IN;
	print OUT "$line[0]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[0]\t",$data[1]/$data[0],"\t",1-$data[-1]/$data[1],"\n";
}
close IL;close OUT;
