#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw/mean median/;
use PerlIO::gzip;
use Getopt::Long;

my $usage="
	Options:
		-in <STR>				The copy ratio file;
		-out <STR>				The output file(average copy ratio);
		-graph <STR>				The output graph (.pdf);
		-cutoff <STR>				The cutoff value for identify aneuploidies [0.75-1.25];
		-R <STR>				The R dir;
		-name <STR>				The file name;
		-h|?					Help!
";
my ($in,$out,$graph,$R,$name,$help);
my $cutoff="0.75-1.25";
GetOptions(
	'in=s' => \$in,
	'out=s' => \$out,
	'graph=s' => \$graph,
	'cutoff=s' => \$cutoff,
	'R=s' => \$R,
	'name=s' => \$name,
	'h|?' => \$help
);

if($help or !$in){die "$usage\n";}
if (-B $in){open IN,"<:gzip","$in";}else{open IN,"<$in";}
open OUT,">$out";
my (%data,@chr_list);
my @cutoff=split /-/,$cutoff;
while(<IN>){
	chomp;
	my @line=split;
	push @{$data{$line[0]}},$line[-1];
}
print OUT "Chr\tCR(Mean)\tCR(Median)\tAneuploidy\tCutoff\n";
open RSH,">$graph.Rsh";
foreach my $chr(sort {$a <=> $b} keys %data){
	print OUT "$chr\t",mean(@{$data{$chr}}),"\t",median(@{$data{$chr}}),"\t",mean(@{$data{$chr}})<$cutoff[0]?"L":mean(@{$data{$chr}})>$cutoff[1]?"G":"N","\t$cutoff\n";
	print RSH "chr$chr=c\(",join (",",@{$data{$chr}}),"\)\n";
	push @chr_list,"chr$chr";
}
print RSH "pdf(\"$graph\",width=16,height=9)\n";
print RSH "par(mar=c(5,5,5,1))\n";
print RSH "boxplot(",join (",",@chr_list),",col=rainbow(24),ylab=\"copy ratio\",xlab=\"chromosome\",",$name?"main=\"$name\",":"","ylim=c(0,3),cex.lab=1.5,cex.main=2,xaxt=\"n\")\n";
print RSH "axis(1,at=c(1:24),labels=c(1:22,\"X\",\"Y\"))\n";
print RSH "dev.off()\n";
print RSH "q()\n";
close IN;close OUT;close RSH;
system "$R CMD BATCH $graph.Rsh $graph.Rout";
system "rm $graph.Rsh $graph.Rout";
