##########################################################################################
### Perform missing value imputation using MICE package
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())
library(mice)
library(abind)

# read in parameters
source(file=commandArgs(trailingOnly=TRUE)[1],echo=T)

# read in metabolite data
metData <- read.table(inFile,header=T,sep="\t",row.names=1,check.names=F)
dim(metData)

sum(colSums(is.na(metData))!=0) # number of mets with NAs

# predictor matrix
corMatrix <- cor(metData,use="pair") # pairwise Pearson correlation, ignore NAs
diag(corMatrix) <- NA

predMatrix <- matrix(0,nrow=dim(metData)[2],ncol=dim(metData)[2])
for (i in 1:dim(corMatrix)[1]) {
	# rank predictor indices by absolute correlation
	rankedCorIndices <- order(abs(corMatrix[i,]),decreasing=T) 
	rankedCors <- corMatrix[i,rankedCorIndices]

	# pick top nPreds predictors, excluding predictors that are too correlated with metabolites (MICE won't take cor > 0.99)
	if (length(which(rankedCors > 0.99)) > 0 ) {
		goodPredIndices <- rankedCorIndices[-which(rankedCors > 0.99)][1:nPreds] 
	} else {
		goodPredIndices <- rankedCorIndices[1:nPreds]
	}
	predMatrix[i,goodPredIndices] <- 1
}

# perform imputation
print("Performing MICE...")
ptm <- proc.time()
imp <- mice(metData, m=nChains, maxit=nIters, method=miceMethod,pred=predMatrix, visit="monotone", print=F)
elapsedT <- proc.time() - ptm
print(paste("MICE TIME:",elapsedT[3]))

# take median of imputed datasets
miceList <- vector(mode="list",length=nChains)
for (iter in 1:nChains) {
	miceList[[iter]] <- complete(imp,iter)
}
medImpData <- apply(abind(miceList,along=3),c(1,2),median,na.rm=T)
dim(medImpData)

sum(colSums(is.na(medImpData))!=0) # number of mets with NAs AFTER imputation (should be few)

# output final imputed data
outTable <- cbind(rownames(metData),signif(medImpData,5))
colnames(outTable)[1] <- "SampleID"
write.table(outTable,outFile,row.name=F,quote=F,sep="\t")

