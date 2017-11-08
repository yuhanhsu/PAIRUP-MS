##########################################################################################
### Filter samples and metabolite signals based on post-QC missingness
### (1) combine QC'ed data from different profiling methods
### (2) remove Internal Standard metabolites
### (3) remove samples with > X (fraction) missingness
### (4) remove metabolite signals with > X missingness
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(plyr) # for rbind.fill()

# read in parameters
source(file=commandArgs(trailingOnly=T)[1], echo=T)

# ----------------------------------------------------------------------------
# (1) combine QC'ed data from different profiling methods

# read in data from each method and store as data frame in a list
dfList <- list()
for (i in 1:length(qcFiles)) {
	dfList[[i]] <- read.table(qcFiles[i],header=T,sep="\t",check.names=F)
	print(dim(dfList[[i]]))
}

# merge data frames based on sample IDs (samples missing in a method woulld have all NAs)
metData <- rbind.fill(dfList)
dim(metData)

# assign row names (metabolite IDs) of combined data frame
rownames(metData) <- unlist(lapply(dfList, row.names))

# ----------------------------------------------------------------------------
# (2) remove Internal Standard metabolites

metData <- metData[-which(rownames(metData) %in% isMets),]
dim(metData)

# transpose met data to get sample x metabolite matrix
sampMetData <- t(metData)
dim(sampMetData) 

# ----------------------------------------------------------------------------
# (3) remove samples with > X missingness

sampMiss <- rowSums(is.na(sampMetData))/dim(sampMetData)[2]
summary(sampMiss)

sampMiss[which(sampMiss > missThreshold)] # samples to remove

if (sum(sampMiss > missThreshold) > 0) {
	sampMetData <- sampMetData[-which(sampMiss > missThreshold),]
}
dim(sampMetData) 

# ----------------------------------------------------------------------------
# (4) remove metabolite signals with > X missingness

metMiss <- colSums(is.na(sampMetData))/dim(sampMetData)[1]
summary(metMiss)

length(metMiss[which(metMiss > missThreshold)]) # number of signals to remove
names(metMiss)[which(metMiss > missThreshold)] # names/IDs of signals to remove

if (sum(metMiss > missThreshold) > 0) {
	sampMetData <- sampMetData[,-which(metMiss > missThreshold)]
}
dim(sampMetData)

# ----------------------------------------------------------------------------
# OUTPUT FINAL DATA

write.table(signif(sampMetData,digits=4),outFile,sep="\t",quote=F)

