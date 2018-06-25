#source("https://bioconductor.org/biocLite.R")
#biocLite("IdeoViz")
#biocLite("IRanges")

library(RColorBrewer)
library(parallel)
library(BiocGenerics)
library(Biobase)
library(stats4)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(rtracklayer)
library(IdeoViz)


PGSoutput = read.table("*.tagsinfo",header = FALSE,sep="\t") ## output file of PGS pipeline.

for (i in 1:nrow(PGSoutput)){
  PGSoutput[i,1]=gsub("23","X",PGSoutput[i,1])
  PGSoutput[i,1]=gsub("24","Y",PGSoutput[i,1])
  PGSoutput[i,1] = paste("chr", PGSoutput[i,1],sep="" ) 

}
ideo_hg19 <- getIdeo("hg19")
chrom_bins <- getBins(paste("chr",c(1:22,"X","Y"), sep=""), ideo_hg19,stepSize=5*100*1000)
# default binning
mean_peak <- avgByBin(data.frame(value=PGSoutput[,7]), PGSoutput[,1:3], chrom_bins)
binned_fullGenome = mean_peak[,2]

pdf(file="IdeoViz2PGS.pdf",pointsize = 10)
plotOnIdeo(chrom=seqlevels(binned_fullGenome),ideo=ideo_hg19,values_GR=binned_fullGenome, 
           value_cols=colnames(mcols(binned_fullGenome)),
           plotType='rect',ablines_y=1,
           col="orange", addScale=F, # hide scale to remove visual clutter
           plot_title="Whole genome view",
           val_range=c(0,2),
           cex.axis=0.6,
           chromName_cex=0.7)
dev.off()
