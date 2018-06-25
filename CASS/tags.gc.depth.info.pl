#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw/mean median/;

use PerlIO::gzip;
use Getopt::Long;

my $usage="
	Options:
		-i/il <STR>				The ext file or list;
		-o <STR>				The output file;
		-w <STR>				The file include of windows[true];
		-h|?					Help!
";
my ($in,$in_list,$out,$help);
my $win="windows.gz";
GetOptions(
	'i=s' => \$in,
	'il=s' => \$in_list,
	'o=s' => \$out,
	'w=s' => \$win,
	'h|?' => \$help
);

if($help or !$out){die "$usage\n";}
if($in_list){
	if(!-d $out){`mkdir -p $out`;}
	open IL,"<$in_list";
	my $out_list=$out."/tags.gc.list";
	open OL,">$out_list";
	while(<IL>){
		chomp;
		my ($name,$ext_path)=(split)[0,-1];
		unless(-f $ext_path){next;}
		my $oo=$out."/$name.gc_cor.gz";
		&tags($ext_path,$oo,$win);
		print OL "$name\t$oo\n";
	}
	close IL;
	close OL;
}
if($in){
	&tags($in,$out,$win);
	#print "Input:$in\nOutput:$out\n$auto_mean\t$genome_mean\t$soap_gc\n";
}
sub tags{
	my ($tin,$tout,$twin)=@_;
	if(-B $tin){open IN,"<:gzip","$tin";}else{open IN,"<$tin";}
	if(-B $twin){open WIN,"<:gzip","$twin";}else{open WIN,"<$twin";}
	open OUT,">:gzip","$tout";
#	open OUT2,">:gzip","$tout.CN.gz";
	my (@ext,@lines,%box,%box_mean);
	my $chr_last=0;
	my $auto_mean=0;my $auto_win_num=0;
#	my $genome_mean=0;my $genome_win_num=0;
	my $reads_num=0;
	while(<WIN>){
		chomp;
		my @win_line=split;
		if(/^#/){next;}
#		next if $win_line[0]>22;
		my $tags=0;
		my $gc=0;
		my $length=0;
		my $bool=1;
		if($chr_last>$win_line[0]){
			$bool=0;
		}
		while($bool){
			my $ext_line=<IN>;
			unless($ext_line){last;}
			chomp($ext_line);
			$reads_num++;
			my @ext_line=split /\s+/,$ext_line;
			$chr_last=$ext_line[0];
			if($ext_line[0]==$win_line[0]){
				if($ext_line[1]<$win_line[2]){
					push @ext,$ext_line;
				}else{
					push @ext,$ext_line;
					last;
				}
			}elsif($ext_line[0]<$win_line[0]){
				undef @ext;
			}else{
				push @ext,$ext_line;
				last;
			}
		}
		while(1){
			last unless $ext[0];
			my ($echr,$elocus)=(split /\s+/,$ext[0])[0,1];
			if($echr<$win_line[0]){
				shift @ext;
			}elsif($echr==$win_line[0]){
				if($elocus<$win_line[1]){shift @ext;}else{last;}
			}else{last;}
		}
		foreach my $e(@ext){
			my ($ext_chr,$ext_locus,$ext_gc,$ext_length)=(split /\s+/,$e)[0,1,2,3];
			if($ext_chr==$win_line[0]){
				if($ext_locus<=$win_line[2] and $ext_locus>=$win_line[1]){
					$tags++;
					$length+=$ext_length;
					$gc+=$ext_gc;
				}
			}
		}
		if($tags==0){print "Null window:$win_line[0]:$win_line[1]-$win_line[2]\n";}
		my $gc_spe=$tags==0?"$win_line[3]":$gc/$length;
		my $win_depth=$length/($win_line[2]-$win_line[1]+1);
		my $lines="$win_line[0]\t$win_line[1]\t$win_line[2]\t$win_line[3]\t$gc_spe\t$win_depth\t$tags";
		push @lines,$lines;
		push @{$box{sprintf("%.2f",$win_line[3])}{sprintf("%.2f",$gc_spe)}},$tags if $win_line[0]<=22;
		if($win_line[0]<=22){
			$auto_win_num++;
			$auto_mean+=$tags;
		}
	}
	my $cor_sum=0;
	$auto_mean=$auto_mean/$auto_win_num;
	foreach(@lines){
		my $gc=sprintf("%.2f",(split)[4]);
		my $ref_gc=sprintf("%.2f",(split)[3]);
		if(!exists $box_mean{$ref_gc}{$gc}){
			if(!exists $box{$ref_gc}{$gc}){
				$box_mean{$ref_gc}{$gc}=$auto_mean;
			}else{
				print $ref_gc,"\t",$gc,"\t",$auto_mean/median(@{$box{$ref_gc}{$gc}}),"\n" unless median(@{$box{$ref_gc}{$gc}})==0;
				$box_mean{$ref_gc}{$gc}=median(@{$box{$ref_gc}{$gc}});
				$box_mean{$ref_gc}{$gc}=$auto_mean if (median(@{$box{$ref_gc}{$gc}})==0 or $auto_mean/median(@{$box{$ref_gc}{$gc}})>10 or $auto_mean/median(@{$box{$ref_gc}{$gc}})<0.1);
			}
		}
		$cor_sum+=((split)[-1]*$auto_mean)/$box_mean{$ref_gc}{$gc} if (split)[0]<=22;
	}
	$cor_sum/=$auto_win_num;
	foreach my $l(@lines){
		print OUT "$l\t",(split /\s+/,$l)[-1]/$auto_mean,"\t",((split /\s+/,$l)[-1]*$auto_mean/$box_mean{sprintf("%.2f",(split /\s+/,$l)[3])}{sprintf("%.2f",(split /\s+/,$l)[4])})/$cor_sum,"\n";
#		print OUT2 "$l\t",(split /\s+/,$l)[-1]/$auto_mean*2,"\t",((split /\s+/,$l)[-1]*$auto_mean/$box_mean{sprintf("%.2f",(split /\s+/,$l)[-3])}{sprintf("%.2f",(split /\s+/,$l)[-2])})/$cor_sum*2,"\n";
	}
	close IN;
	close OUT;
#	close OUT2;
	close WIN;
#	$genome_mean=$genome_mean/$genome_win_num;
}
