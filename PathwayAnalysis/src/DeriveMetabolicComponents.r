##########################################################################################
### (1) Calculate pairwise Spearman correlation matrix for metabolite data
### (2) Perform PCA on correlation amtrix to derive "Metabolic Components" (MC)
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# read in parameters from input argument file name
source(commandArgs(trailingOnly=T)[1],echo=T)

# ----------------------------------------------------------------------------------------
# read in metabolite data
metData <- read.table(metDataFile,sep="\t",header=T,row.names=1,check.names=F)
dim(metData)

# calculate correlation matrix
print("Calculating correlation matrix:")
corMatrix <- cor(metData, method="spearman")
dim(corMatrix)

# output correlatrion matrix (optional)
if (length(corMatrixFile) > 0) {
	write.table(signif(corMatrix,digits=4),file=gzfile(corMatrixFile),sep="\t",quote=F)
}

# set any missing values in correlation matrix to 0
corMatrix[is.na(corMatrix)] <- 0

# perform principal component analysis
print("Performing PCA:")
pcaResults <- prcomp(corMatrix)

# output stats for 1st nSamples MCs
# (standard deviation, proportion of variance, cumulative proportion of variance)
SD <- pcaResults$sdev[1:nSamples]
PropVar <- SD^2 / sum(SD^2)
CumProp <- cumsum (PropVar)
outDf <- data.frame(MC=1:nSamples, StandardDeviation=SD, ProportionVariance=PropVar, CumulativeProportion=CumProp)
write.table(signif(outDf,digits=4),file=mcStatsFile,sep="\t",quote=F,row.names=F)

# output signal x MC (PC) score matrix for 1st nSamples MCs
outTable <- signif(pcaResults$x[,1:nSamples],digits=4)
colnames(outTable) <- paste("MC",1:nSamples,sep="")
write.table(outTable,file=mcScoreMatrixFile,sep="\t",quote=F)

print("SCRIPT COMPLETED")
