####by Heather E. Wheeler 20160421####
###INPUT: PrediXcan dosage files
###OUTPUT: MACH dosage and mlinfo files for gcta
date <- Sys.Date()
"%&%" = function(a,b) paste(a,b,sep="")
args <- commandArgs(trailingOnly=T)

###############################################
### Directories & Variables

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
exp.dir <- my.dir %&% "GEUVADIS/"
gt.dir <- exp.dir %&% "geu_dosage_1KGP_MAF0.01/"

################################################
### Functions & Libraries
library(dplyr)
library(data.table)
################################################

samples <- read.table(gt.dir %&% "samples.txt")
sampleinfo <- dplyr::select(samples,V1) %>% mutate(V1=paste(V1,"->",V1,sep=""),V2="MLINFO")

i <- args[1]

gtfile <- "gunzip -c " %&% gt.dir %&% "chr" %&% i %&% ".dosage.gz"
gt <- fread(gtfile)
gtinfo <- gt[,1:6,with=FALSE]
gtmat <- as.matrix(gt[,7:dim(gt)[2],with=FALSE])

colnames(gtinfo) <- c('chr','SNP','pos','freq2','Al1','Al2')
gtinfo <- mutate(gtinfo, Freq1=round((1-freq2),3), MAF=ifelse(freq2>0.5,round(1-freq2,3),round(freq2,3)))
mlinfo <- dplyr::select(gtinfo, SNP, Al1, Al2, Freq1, MAF) %>% mutate(Quality=1, Rsq=1)

write.table(mlinfo, gt.dir %&% "chr" %&% i %&% ".mlinfo",quote=F,row.names=F)

#transpose for mach format
tgtmat <- t(gtmat)
#add first 2 cols of mach mldose format
mldose <- cbind(sampleinfo, tgtmat)
write.table(mldose, gt.dir %&% "chr" %&% i %&% ".mldose", quote=F, row.names = F, col.names = F)
  
#compress files
runGZIP <- "gzip " %&% gt.dir %&% "chr" %&% i %&% ".ml*"
system(runGZIP)

