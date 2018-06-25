#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw/mean median/;
use PerlIO::gzip;
use Getopt::Long;

my $usage="
	Options:
		-i <STR>		The copy ratio file or list(split by ',');
		-out <STR>		The output folder;
		-name <STR>		The name of data;
		-R <STR>		The R dir[/usr/bin/R];
		-h|?			Help!
";
my ($in,$out,$name,$help);
my $R="/usr/bin/R";
GetOptions(
	'i=s' => \$in,
	'out=s' => \$out,
	'name=s' => \$name,
	'R=s' => \$R,
	'h|?' => \$help
);

if($help or !$out){die "$usage\n";}
gc_CR();
GC_depth();
sub gc_CR{
	my @in=split /,/,$in;
	die "Wrong input!\n" if @in>2;
	my @name=split /,/,$name;
	open SH,">$out/GC_CR.Rsh";
	print SH "pdf(\"$out/$name[-1]_GC_CR.pdf\",width=16,height=9)\n";
	my $line_num=0;
	my @col=("rgb(0,0,1,0.8)","rgb(1,0,0,0.3)");
	my $samples="\"".join("\",\"",@name)."\"";
	foreach my $i(@in){
		$line_num++;
		print SH $line_num==1?"par(mar=c(5,5,1,1))\n":"par(new=TRUE)\n";
		if(-B $i){open IN,"<:gzip","$i";}else{open IN,"<$i";}
		my (@gc,@ratio);
		while(<IN>){
			chomp;
			my @line=split;
			push @gc,$line[4];
			push @ratio,$line[7];
		}
		close IN;
		my $gc_column=join ",",@gc;
		my $ratio_column=join ",",@ratio;
		print SH "GC=c($gc_column)\n";
		print SH "RT=c($ratio_column)\n";
		print SH "GC_RT=cbind(GC,RT)\n";
		print SH "plot(GC_RT,type=\"p\",pch=20,col=$col[$line_num-1],cex=0.5,xlim=c(0.2,0.8),ylim=c(0,3),",$line_num==1?"xlab=\"GC\",ylab=\"Copy Ratio\",":"xlab=\"\",ylab=\"\",","cex.lab=1.5,xaxt=\"n\")\n";
		print SH "axis(1,at=seq(0,1,0.1))\n" if $line_num==1;
	}
	my $legend_col=$line_num==2?(join ",",@col):$col[0];
	print SH "legend(\"topright\",c($samples),col=c($legend_col),pch=20,cex=1.5,bty=\"n\")\n";
	print SH "dev.off()\nq()\n";
	close SH;
	system "$R CMD BATCH $out/GC_CR.Rsh $out/GC_CR.Rout";
}
sub GC_depth{
	my @in=split /,/,$in;
	die "Wrong input!\n" if @in>2;
	my @name=split /,/,$name;
	open SH,">$out/GC_depth.Rsh";
	my $line_num=0;
	foreach my $i(@in){
		if(-B $i){open IN,"<:gzip","$i";}else{open IN,"<$i";}
		$line_num++;
		my (@chr_locus,@chr_gc,@chr_depth);
		my $last_chr=1;
		while(<IN>){
			chomp;
			my @line=split;
			if($line[0]!=$last_chr){
				print SH "pdf(\"$out/$name[$line_num-1]_chr$last_chr.GC_depth.pdf\",width=16,height=8)\n";
				print SH "par(mar=c(5,5,1,5))\n";
				print SH "GC$last_chr=c(",join(",",@chr_gc),")\n";
				print SH "Locus$last_chr=c(",join(",",@chr_locus),")\n";
				print SH "Depth$last_chr=c(",join(",",@chr_depth),")\n";
				print SH "GC_Locus$last_chr=cbind(Locus$last_chr,GC$last_chr)\n";
				print SH "Depth_Locus$last_chr=cbind(Locus$last_chr,Depth$last_chr)\n";
				print SH "plot(Depth_Locus$last_chr,type=\"h\",lwd=1,col=rgb(0,0,1),cex=1,ylim=c(0,10),xlab=\"ChromosomeLocus\",ylab=\"Depth\",cex.lab=1.5)\n";
				print SH "par(new=TRUE)\n";
				print SH "plot(GC_Locus$last_chr,type=\"p\",pch=20,col=rgb(1,0,0),cex=0.5,ylim=c(0,1),xlab=\"\",ylab=\"\",cex.lab=1.5,yaxt=\"n\")\n";
				print SH "axis(4,at=seq(0,1,0.1))\n";
				print SH "mtext(\"GC\",side=4,cex=1.5,line=3)\n";
				print SH "legend(\"topright\",c(\"Depth\",\"GC\"),col=c(rgb(0,0,1),rgb(1,0,0)),pch=20,cex=1,bty=\"n\")\n";
				print SH "dev.off()\n";
				undef @chr_locus;
				undef @chr_gc;
				undef @chr_depth;
				push @chr_gc,$line[4];
				push @chr_depth,$line[7];
				push @chr_locus,($line[1]+$line[2])/2;
				$last_chr=$line[0];
				next;
			}
			push @chr_gc,$line[4];
			push @chr_depth,$line[7];
			push @chr_locus,($line[1]+$line[2])/2;
			$last_chr=$line[0];
		}
		print SH "pdf(\"$out/$name[$line_num-1]_chr$last_chr.GC_depth.pdf\",width=16,height=8)\n";
		print SH "par(mar=c(5,5,1,5))\n";
		print SH "GC$last_chr=c(",join(",",@chr_gc),")\n";
		print SH "Locus$last_chr=c(",join(",",@chr_locus),")\n";
		print SH "Depth$last_chr=c(",join(",",@chr_depth),")\n";
		print SH "GC_Locus$last_chr=cbind(Locus$last_chr,GC$last_chr)\n";
		print SH "Depth_Locus$last_chr=cbind(Locus$last_chr,Depth$last_chr)\n";
		print SH "plot(Depth_Locus$last_chr,type=\"h\",lwd=1,col=rgb(0,0,1),cex=1,ylim=c(0,10),xlab=\"ChromosomeLocus\",ylab=\"Depth\",cex.lab=1.5)\n";
		print SH "par(new=TRUE)\n";
		print SH "plot(GC_Locus$last_chr,type=\"p\",pch=20,col=rgb(1,0,0),cex=0.5,ylim=c(0,1),xlab=\"\",ylab=\"\",cex.lab=1.5,yaxt=\"n\")\n";
		print SH "axis(4,at=seq(0,1,0.1))\n";
		print SH "mtext(\"GC\",side=4,cex=1.5,line=3)\n";
		print SH "legend(\"topright\",c(\"Depth\",\"GC\"),col=c(rgb(0,0,1),rgb(1,0,0)),pch=20,cex=1,bty=\"n\")\n";
		print SH "dev.off()\n";
	}
	print SH "q()\n";
	close SH;
	system "$R CMD BATCH $out/GC_depth.Rsh $out/GC_depth.Rout";
}
