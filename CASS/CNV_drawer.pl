#!/usr/bin/perl -w
use strict;
use SVG;
use Getopt::Long;
use File::Basename;
use Cwd qw/abs_path/;

my $usage=<<'USAGE';
Usage: 

	perl CNV_Drawer.pl [Options]

	-i FILE    CNV result
	-o DIR     output folder
	-z FILE    the correction ratio
	-r EXE     the R dir
	-ylim STR  the ylim during spots plotting [0,2];
	-name STR  the file name;
	-h HELP    help information


USAGE

my ($in,$out,$ratio,$R,$name,$help);
my $ylim="0,2";
#my @chr=(247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,134452384,132349534,114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954);
my @chr=(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566);
my @centromere=(123035434,93826171,92004854,51160117,47905641,60330166,59554331,45338887,48867679,40754935,53144205,36356694,17500000,17500000,18500000,36835801,23763006,16960898,26181782,27869569,12788129,14500000,60132012,11604553);
GetOptions(
	'i=s'=>\$in,
	'o=s'=>\$out,
	'z=s'=>\$ratio,
	'r=s'=>\$R,
	'ylim=s' =>\$ylim,
	'name=s' => \$name,
	'h|?'=>\$help,
);

die "$usage\n" if($help or !$out);

$in=abs_path($in);
$out=abs_path($out);
&spot($in,$out,$ratio);
sub spot{
    my ($in,$dir,$ratio)=@_;
    if(!-d $dir){`mkdir -p $dir`;}
    open SH,">$dir/draw.Rsh";
    if (-B $ratio){open ZS,"<:gzip","$ratio";}else{open ZS,"<$ratio";}
    my $chr=0;
    my (%data,%locus,@genome_data,@genome_locus);
    while(<ZS>){
        chomp;
        next unless /^\d+/;
        my @line=split;
        push @{$data{$line[0]}},$line[-1];
        push @{$locus{$line[0]}},($line[2]+$line[1])/2;
		my $length=0;
		for(1..($line[0]-1)){$length+=$chr[$_-1];}
		push @genome_data,$line[-1];
		push @genome_locus,$length+($line[2]+$line[1])/2;
        $chr=$line[0];
    }
    close ZS;
    foreach my $c(sort {$a <=> $b} keys %data){
        my $data=join ",",@{$data{$c}};
        my $locus=join ",",@{$locus{$c}};
        print SH "locus$c=c($locus)\n";
        print SH "zsocre$c=c($data)\n";
        print SH "chr$c=cbind(locus$c,zsocre$c)\n";
    }
	my $genome_data=join ",",@genome_data;
	my $genome_locus=join ",",@genome_locus;
	print SH "locus=c($genome_locus)\n";
	print SH "zscore=c($genome_data)\n";
	print SH "genome=cbind(locus,zscore)\n";
    open IN,"<$in";
	my $last_line="";
    my $chr_num=1;
	my $bool=1;
	my $genome_seq="";
    for(1..$chr){
        my $dchr=$_;
        my $png=$name?$dir."/$name.chr$dchr.pdf":$dir."/chr$dchr.pdf";
		my $seg_length=0;
		for(1..($dchr-1)){$seg_length+=$chr[$_-1];}
		print SH "pdf(\"$png\",width=16,height=12)\n";
        print SH "par(mar=c(5,5,5,1))\n";
        print SH "plot(chr$dchr,type=\"p\",pch=20,,col=\"blue\",cex=0.5,ylim=c($ylim),xlab=\"Locus\",ylab=\"copy ratio\",",$name?"main=\"$name chr$dchr\"":"main=\"chr$dchr\"",",cex.lab=1.5,cex.main=2)\n";
	    if($last_line){
			my @last_line=split /\s+/,$last_line;
			my $last=$last_line[1].",".$last_line[3].",".$last_line[2].",".$last_line[3];
			my $genome_last=($seg_length+$last_line[1]).",".$last_line[3].",".($seg_length+$last_line[2]).",".$last_line[3];
			if($last_line[0]==$dchr){print SH "segments($last,lwd=3,col=\"red\")\n";$genome_seq.="segments($genome_last,lwd=3,col=\"red\")\n";}else{$bool=0;}
		}
    	while($bool){
            my $dna_line=<IN>;
           	last unless $dna_line;
	        next unless $dna_line=~/^\d+/;
    	    my @dna_line=split /\s+/,$dna_line;
            if($dna_line[0]==$dchr){
            	my $seg=$dna_line[1].",".$dna_line[3].",".$dna_line[2].",".$dna_line[3];
				my $genome_seg=($seg_length+$dna_line[1]).",".$dna_line[3].",".($seg_length+$dna_line[2]).",".$dna_line[3];
                print SH "segments($seg,lwd=3,col=\"red\")\n";
				$genome_seq.="segments($genome_seg,lwd=3,col=\"red\")\n";
	        }else{
				$last_line=$dna_line;
        	    last;
            }
		}
		print SH "dev.off()\n";
		$chr_num=$dchr;
		#print SH "legend(\"topright\",c(\"copy ratio\",\"Segmentation\"),col=c(\"blue\",\"red\"),lwd=2)\n";
	}
	my $genome_png=$name?$dir."/$name.genome.pdf":$dir."/genome.pdf";
    print SH "pdf(\"$genome_png\",width=16,height=9)\n";
    print SH "par(mar=c(5,5,5,1))\n";
    print SH "plot(genome,type=\"p\",pch=20,,col=\"blue\",cex=0.5,ylim=c($ylim),xlab=\"Locus\",ylab=\"copy ratio\",",$name?"main=\"$name genome\"":"main=\"genome\"",",cex.lab=1.5,cex.main=2,xaxt=\"n\")\n";
	print SH "$genome_seq";
	my $chr_length=0;
	my $chr_locus="";
	my @tick;
	for(1..24){$chr_length+=$chr[$_-1];$chr_locus.=",".$chr_length;push @tick,(($chr_length+$chr_length-$chr[$_-1])/2);}
	my $tick=join ",",@tick;
	$chr_locus="0".$chr_locus;
	print SH "axis(1,at=c($chr_locus),labels=FALSE)\n";
	print SH "axis(1,at=c($tick),labels=c(1:22,\"X\",\"Y\"),tick=FALSE)\n";
	print SH "dev.off()\n";
    print SH "q()\n";
    close SH;close IN;
    system "$R CMD BATCH $dir/draw.Rsh $dir/draw.Rout";
#	system "rm $dir/draw.Rsh $dir/draw.Rout";
}
