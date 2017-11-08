##########################################################################################
### Adjust and transform metaoblite data after missing value imputation
### (1) perform median imputation for the few NAs that were not imputed by MICE
### (2) adjust for covariates (if specified in param file)
### (3) perform rank-based inverse normal transformation
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# read in parameters
source(file=commandArgs(trailingOnly=TRUE)[1],echo=T)

# read in metabolite data
metData <- read.table(metFile,sep="\t",header=T,row.names=1,check.names=F)
dim(metData)

# --------------------------------------------------------------------------------------
# (1): fill in missing data (should be few after MICE) with median of metabolite value
sum(is.na(metData))

metMeds = apply(metData,2,median,na.rm=TRUE)
for (i in 1:length(metMeds)) {
	metData[which(is.na(metData[,i])),i] = metMeds[i]
}

sum(is.na(metData))

# -------------------------------------------------------------------------------------
# (2): adjust for covariates (optional)
# (3): inverse normal transformation using the residuals

# get covariate data and remove samples with missing covariates
if (length(covariates) > 0) { # if any covariates were specified in param file
	
	# read in phenotype data
	ptypeData <- read.table(ptypeFile,header=T,sep="\t",row.names=1)
	
	# sort phenotype table to have same sample order as metData
	ptypeTable <- ptypeData[match(rownames(metData),rownames(ptypeData)),]

	# remove samples with missing covariates info
	badCovSamples <- rep(FALSE,dim(ptypeTable)[1])
	for (i in covariates) {
		badCovSamples <- badCovSamples | is.na(ptypeTable[,i])
	}

	if (sum(badCovSamples) > 0) { # need to remove samples
		covTable <- ptypeTable[-which(badCovSamples),covariates]
		metData <- metData[-which(badCovSamples),]
	} else { # if no samples to remove
		covTable <- ptypeTable[,covariates]
	}
	dim(covTable)
	dim(metData)
}

outTable <- NULL
for (i in 1:dim(metData)[2]) {
	m = metData[,i]
	
	if (length(covariates) > 0) { # if adjusting for covariates
		model = lm(m ~ ., data=covTable, na.action=na.exclude)
		res = residuals(model)
	} else { # no adjustment
		res = m
	}

	# rank-based inverse normal transformation
	transM = qnorm((rank(res, na.last="keep") - 0.5) / sum(!is.na(res)))

	# bind new metabolite column to output table
	outTable <- cbind(outTable,signif(transM,digits=4))
}

colnames(outTable) <- colnames(metData)
rownames(outTable) <- rownames(metData)

# output adjusted and transformed metabolite data
write.table(outTable, file=outFile, quote=F, sep="\t")

