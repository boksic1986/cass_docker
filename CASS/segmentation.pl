#!/usr/bin/perl -w

## For: CNV calling for PGS based on NGS;
## Author: Chunlei Zhang;
## Email: yxtkttkl@gmail.com;
## Version: v14.10

use strict;
use Statistics::Distributions;
use Statistics::Sequences::Runs;
use Statistics::Data;
use Getopt::Long;
use Statistics::Basic qw/mean median/;
use List::MoreUtils qw/uniq/;
use List::Util qw/sum min max/;

my $usage=<<'USAGE';
Usage: 

    perl Segmentation.pl [Options]

    -i FILE    	ratio/z-score file
    -o FILE    	output file
    -n FILE    	path of N region in the chromosomes [hg18]
    -w INT     	number of windows in each initial test [100]
    -p FLOAT   	threshold of the initial test [1]
    -d INT     	number of pionts to be retained in the initial test [10000]
    -e FLOAT   	finish threshold of the secend test [1e-25]
    -fa FLOAT	the lower cutoff value for identify deletions [0.75]
    -fb FLOAT	the upper cutoff value for identify duplications [1.25]
    -gender STR	the gender info [false]
	-band STR	the cyto band for annotation
    -h HELP    help information

Example

    perl Segmentation.pl -i ./Sample1.ratio -o Sample1.cnv
    perl Segmentation.pl -i ./Sample2.ratio -p 1e-10 -o Sample2.cnv
    perl Segmentation.pl -i ./Sample3.ratio -p 1e-3 -d 5000 -o Sample3.cnv
    perl Segmentation.pl -i ./Sample4.ratio -p 1e-3 -d 5000 -e 1e-30 -o Sample4.cnv
USAGE

my $N_region = "/ifs1/ST_REHEAL/USER/zhangchl/soft/cnv/nips/FCAPS/n_region/hg18.cen_ex1M.N";
my ($in,$out,$gender,$cytoband,$help);
my ($window,$initial_p,$initial_d,$end_p,$percent,$fa,$fb) = (100,1,10000,"1e-25",0.015,0.75,1.25);

GetOptions(
    'in=s' => \$in,
	'out=s' => \$out,
	'n_region=s' => \$N_region,
	'window=i' => \$window,
	'p=s' => \$initial_p,
	'd=i' => \$initial_d,
	'e=s' => \$end_p,
	'fa=s' => \$fa,
	'fb=s' => \$fb,
	'gender=s' => \$gender,
	'band=s' => \$cytoband,
	'help|?' => \$help,
);

die "$usage\n" if ($help || !$in || !$out);

#-------------------------cytoband annotation (2015.06.12)---------------------------------
my %cytoband=readband($cytoband) if $cytoband;
#-------------------------cytoband annotation (2015.06.12)---------------------------------
if(-B $in){open IN,"<:gzip","$in";}else{open IN,"<$in";}
open NR,"<$N_region";
open OUT,">$out";
print OUT "chr\tstart\tend\tcopy ratio/number\tsize\t",$cytoband?"cytoband\t":"","p-value\n";

my ($pre_chr,$pre_pos,$pre_score,$chr_length,%Chromosome,%N_region,%Score,%BreakPiont,$breakpiont);
my $zed = Statistics::Zed -> new();

print "read N_region ...\n";
while(<NR>){
	chomp;
	my ($chr,$start,$end) = split;
	$chr=~s/chr//g;
	$chr=($chr=~/X/)?23:($chr=~/Y/)?24:$chr;
	$N_region{$chr}{$start} = $end;
}
close NR;

print "read ratio of case ...\n";
while(<IN>){
	chomp;
	next unless /^\d/;
	my ($chr,$start,$end,$score) = (split)[0,1,2,-1];
	if($chr>22){$score=($gender and $gender eq "M")?($score*2):$score;}
	my $bool = 0;
	foreach my $s(keys %{$N_region{$chr}}){
		my $e = $N_region{$chr}{$s};
		$bool = 1 if ($s < $end && $e > $start);
	}
	next if $bool;
	$Chromosome{$chr}{"s"} = $chr * 1e9 + $start unless exists $Chromosome{$chr}{"s"};
	$Chromosome{$chr}{"e"} = $chr * 1e9 + $end;
	my $pos = $chr * 1e9 + int(($start + $end)/2);
	$Score{$pos} = $score;
}
close IN;

print "Initial P ...\n";
my @score = sort {$a <=> $b} keys %Score;
foreach my $index(0..$#score){
	my @left = map {$Score{$_}} @score[$index-$window..$index];
	my @right = map {$Score{$_}} @score[map {if($_ > $#score){$_ -= scalar @score}else{$_}} ($index..$index+$window)];
#	my $test1=join "\t",(sort {$a <=> $b} @left);my $test2=join "\t",(sort {$a <=> $b} @right);die "$test1\n****************\n$test2\n";
	my $p_value = Runs_Test(\@left,\@right);
	if($p_value < $initial_p){
		$BreakPiont{$score[$index]} = $p_value;
		$breakpiont ++;
	}
}

if($breakpiont > $initial_d){
	my $rank_d = (sort {$a<=> $b} values %BreakPiont)[$initial_d-1];
	map {delete $BreakPiont{$_} if $BreakPiont{$_} > $rank_d} keys %BreakPiont;
	$breakpiont = scalar keys %BreakPiont;
}

my @breakpiont = sort {$a <=> $b} keys %BreakPiont;

foreach my $index(0..$#breakpiont){
	my $pos = $breakpiont[$index];
	my ($left,$right) = @breakpiont[$index-1,$index+1];
	$right ||= $breakpiont[0];
	my @left = $left < $pos ? map {$Score{$_}} grep {$_ >= $left && $_ <= $pos} @score : map {$Score{$_}} grep {$_ >= $left || $_ <= $pos} @score;
	my @right = $pos < $right ? map {$Score{$_}} grep {$_ >= $pos && $_ <= $right} @score : map {$Score{$_}} grep {$_ >= $pos || $_ <= $right} @score;
	$BreakPiont{$pos} = Runs_Test(\@left,\@right);
}

print "Segmentation from $breakpiont ...\n";
my ($mutation,$cycle) = (0,0);

while(1){
	my ($key,$value) = (0,0,0);
	foreach my $hash_key(keys %BreakPiont){
		if($BreakPiont{$hash_key} > $value){
			$value = $BreakPiont{$hash_key};
			$key = $hash_key;
			last if $value == 1;
		}
	}
	$cycle ++;
	
	my $left = $breakpiont - $cycle;
	last if ($value <= $end_p || $left <= 24);

	print "cycle $cycle : $left $value ...\n" if $cycle % 100==1;

	my @pos = sort {$a <=> $b} keys %BreakPiont;

	if(scalar @pos > 1){
		my ($index) = grep {$pos[$_] == $key} (0..$#pos);
		my ($left_left,$left) = @pos[$index-2,$index-1];
		my ($right,$right_right) = @pos[$index+1,$index+2];
		if(!defined $right){
			$right = $pos[0];
			$right_right = $pos[1];
		}elsif(!defined $right_right){
			$right_right = $pos[0];
		}

		my @left_left = $left_left < $left ? map {$Score{$_}} grep {$_ >= $left_left && $_ <= $left} @score : map {$Score{$_}} grep {$_ >= $left_left || $_ <= $left} @score;
		my @middle = $left < $right ? map {$Score{$_}} grep {$_ >= $left && $_ <= $right} @score : map {$Score{$_}} grep {$_ >= $left || $_ <= $right} @score;
		my @right_right = $right < $right_right ? map {$Score{$_}} grep {$_ >= $right && $_ <= $right_right} @score : map {$Score{$_}} grep {$_ >= $right || $_ <= $right_right} @score;
	
		$BreakPiont{$left} = Runs_Test(\@left_left,\@middle);
		$BreakPiont{$right} = Runs_Test(\@middle,\@right_right);
	}

	delete $BreakPiont{$key};
}

print "output result ...\n";
foreach my $chr(sort {$a <=> $b} uniq(map {int($_ / 1e9)} keys %Score)){
#	if($gender ne "M" and $chr==24){next;}
    if(my @pos = sort {$a <=> $b} grep {$_ >= $chr * 1e9 && $_ <= ($chr + 1) * 1e9} keys %BreakPiont){
        foreach my $index(0..$#pos){
            my $left = $index ? $pos[$index-1] : $Chromosome{$chr}{"s"};
            my $right = $pos[$index+1] || $Chromosome{$chr}{"e"};
			my @cutn=cutn($chr,$left-$chr * 1e9,$pos[$index]- $chr * 1e9,\%N_region);
			
			for(my $i=0;$i<=$#cutn;$i+=3){
	            my @left_score = map {$Score{$_}} grep {$_ >= $cutn[$i+1]+$chr * 1e9 && $_ <= $cutn[$i+2]+ $chr * 1e9} keys %Score;
        	    next unless @left_score;
#            	my ($left_ratio) = Mean_SD(@left_score);
            	my ($left_ratio) = median(@left_score);
    	        my $left_p_value = pvalue(@left_score);

        	    my $left_support = scalar @left_score;
	            print OUT "$chr\t$cutn[$i+1]\t",$cutn[$i+2],"\t",sprintf("%.3f",$left_ratio);
				print OUT "\t",$cutn[$i+2]-$cutn[$i+1]>=1000000?((sprintf("%.2f",($cutn[$i+2]-$cutn[$i+1])/1000000))."Mb"):$cutn[$i+2]-$cutn[$i+1]>=1000?((sprintf("%.2f",($cutn[$i+2]-$cutn[$i+1])/1000))."Kb"):($cutn[$i+2]-$cutn[$i+1]);
				print OUT "\t",$left_ratio>$fb?"G":$left_ratio<$fa?"L":"N";
				print OUT "\t",bandanno($chr,$cutn[$i+1],$cutn[$i+2],\%cytoband) if $cytoband;
				print OUT "\t",sprintf("%.2e",$left_p_value),"\n";
			}

            if($right == $Chromosome{$chr}{"e"}){
				@cutn=cutn($chr,$pos[$index]-$chr * 1e9,$right-$chr * 1e9,\%N_region);
				for(my $j=0;$j<=$#cutn;$j+=3){
	                my @right_score = map {$Score{$_}} grep {$_ >= $cutn[$j+1]+$chr * 1e9 && $_ <= $cutn[$j+2]+$chr * 1e9} keys %Score;
        	        next unless @right_score;
	
#    	            my ($right_ratio) = Mean_SD(@right_score);
    	            my ($right_ratio) = median(@right_score);
					my $right_p_value = pvalue(@right_score);
	
    	            my $right_support = scalar @right_score;
            	    print OUT "$chr\t",$cutn[$j+1],"\t$cutn[$j+2]\t",sprintf("%.3f",$right_ratio);
					print OUT "\t",$cutn[$j+2]-$cutn[$j+1]>=1000000?((sprintf("%.2f",($cutn[$j+2]-$cutn[$j+1])/1000000))."Mb"):$cutn[$j+2]-$cutn[$j+1]>=1000?((sprintf("%.2f",($cutn[$j+2]-$cutn[$j+1])/1000))."Kb"):($cutn[$j+2]-$cutn[$j+1]);
					print OUT "\t",$right_ratio>$fb?"G":$right_ratio<$fa?"L":"N";
					print OUT "\t",bandanno($chr,$cutn[$j+1],$cutn[$j+2],\%cytoband) if $cytoband;
					print OUT "\t",sprintf("%.2e",$right_p_value),"\n";
				}
            }
        }
    }else{
		my ($left,$right) = ($Chromosome{$chr}{"s"},$Chromosome{$chr}{"e"});
        my @middle_score = map {$Score{$_}} grep {$_ >= $left && $_ <= $right} keys %Score;
		my @middle_break = grep {$_ >= $left && $_ <= $right} keys %Score;

#        my ($middle_ratio) = Mean_SD(@middle_score);
        my ($middle_ratio) = median(@middle_score);
		my $middle_p_value = pvalue(@middle_score);

        my $middle_support = scalar @middle_score;
        $left -= $chr * 1e9;
        $right -= $chr * 1e9;
        print OUT "$chr\t$left\t$right\t",sprintf("%.3f",$middle_ratio);
		print OUT "\t",$right-$left>=1000000?((sprintf("%.2f",($right-$left)/1000000))."Mb"):$right-$left>=1000?((sprintf("%.2f",($right-$left)/1000))."Kb"):($right-$left);
		print OUT "\t",$middle_ratio>$fb?"G":$middle_ratio<$fa?"L":"N";
		print OUT "\t",bandanno($chr,$left,$right,\%cytoband) if $cytoband;
		print OUT "\t",sprintf("%.2e",$middle_p_value),"\n";
    }
}
close OUT;

sub SD{
    return (sum(map {$_ ** 2} @_) / (scalar @_)) ** 0.5;
}
sub pvalue{
    my (@data)=@_;
    my $pvalue=1;
    if(@data>1){
		my ($mean,$sd)=Mean_SD(values %Score);
        my $sum=0;
		foreach (@data){$sum+=(($_-$mean)/$sd)**2;}
        $pvalue=Statistics::Distributions::chisqrprob($#data,$sum);
    }
    $pvalue=sprintf("%.2e",$pvalue);
    return $pvalue;
}
sub cutn{
    my ($chr,$start,$end,$nr)=@_;
    my %nr=%$nr;
    my ($bool1,$bool2)=(1,0);
    my @cut;
    foreach my $s(sort {$a <=> $b} keys %{$nr{$chr}}){
        my $e=$nr{$chr}{$s};
		next if $e-$s+1<200000;
        if($s<=$end and $e>=$start){
            $bool1=0;
            if($s>$start){
                push @cut,$chr;
                push @cut,$start;
                push @cut,($s-1);
                if($end>$e){$bool2=1;$start=$e+1;}else{$bool2=0;}
            }else{
                if($end>$e){$bool2=1;$start=$e+1;}else{$bool2=0;}
            }
        }
    }
    if($bool1 or $bool2){
        push @cut,$chr;
        push @cut,$start;
        push @cut,$end;
    }
    return @cut;
}
sub Open{
	my ($file,$symbol)=@_;
	$symbol ||= "<";
	my $handle;
	if($_[0] =~ /\.gz$/ && ($symbol eq ">" || -B $_[0])){
		open $handle,"$symbol:gzip",$_[0];
	}else{
		open $handle,$symbol,$_[0];
	}
	return $handle;
}
sub Runs_Test{
    my ($array1,$array2)=@_;
    my @array1=@$array1;
    my @array2=@$array2;
    my %data;
    foreach(0..$#array1){
        while(1){
            last unless exists $data{$array1[$_]};
            $array1[$_]+=(rand(1)-0.5)/1e12;
        }
        $data{$array1[$_]}=1;
    }
    foreach(0..$#array2){
        while(1){
            last unless exists $data{$array2[$_]};
            $array2[$_]+=(rand(1)-0.5)/1e12;
        }
        $data{$array2[$_]}=0;
    }
    my @runs;
    foreach(sort {$a <=> $b} keys %data){
        push @runs,$data{$_};
    }
    my $runs = Statistics::Sequences::Runs->new();
    $runs->load(@runs);
    $runs->test();
#   return $runs->{'p_value'};
    return $runs->p_value;
}
sub Median_Test{
	my @group1=@{$_[0]};
	my @group2=@{$_[1]};
	my @group=sort(@group1,@group2);
	my $median=median(@group);
	my $n11=scalar (grep {$_ > $median} @group1);
	my $n1p=scalar (grep {$_ > $median} @group);
	my $np1=scalar @group1;
	my $npp=scalar @group;
	my $left_value = calculateLeftStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	my $right_value = calculateRightStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	return (min($left_value,$right_value)*2 > 1 ? 1 : min($left_value,$right_value)*2);
}
sub Mean_SD{
	die "The array is empty!\n" unless @_;
	@_ = @{$_[0]} if (@_ == 1 && $_[0] =~ /^ARRAY/);
	if(@_ == 1){
		return ($_[0],0);
	}else{
		my $mean = sum(@_)/@_;
		my $sd = (sum((map {($_-$mean)**2} @_))/(@_-1))**0.5;
		return ($mean,$sd);
	}
}
#-------------------------cytoband annotation (2015.06.12)---------------------------------
sub readband{
	if(-B $cytoband){open BAND,"<:gzip","$cytoband";}else{open BAND,"<$cytoband";}
	my %hash;
	while(<BAND>){
		chomp;
		my @line=split;
		$line[0]=~s/chr//ig;
		$line[0]=$line[0]=~/X/?23:$line[0]=~/Y/?24:$line[0]=~/M/?25:$line[0];
		$hash{$line[0]}{$line[1]}{$line[2]}=$line[3];
	}
	close BAND;
	return %hash;
}
sub bandanno{
	my ($bchr,$bstart,$bend,$hash)=@_;
	my %bhash=%$hash;
	my ($startband,$endband);
	foreach my $s(sort {$a <=> $b} keys %{$bhash{$bchr}}){
		foreach my $e(sort {$a <=> $b} keys %{$bhash{$bchr}{$s}}){
			$startband=$bhash{$bchr}{$s}{$e} if ($bstart>=$s and $bstart<=$e);
			$endband=$bhash{$bchr}{$s}{$e} if ($bend>=$s and $bend<=$e);
		}
	}
	if($startband eq $endband){return "$startband";}else{return "$startband-$endband";}
}
#-------------------------cytoband annotation (2015.06.12)---------------------------------
