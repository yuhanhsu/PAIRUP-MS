##########################################################################################
### Impute signals across datasets using shared known metabolites as predictors
### (shared knowns were imputed using leave-one-out models)
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# TIME PROCESS
ptm <- proc.time()

# read in parameters from input argument file name
source(commandArgs(trailingOnly=T)[1],echo=T)


# ----------------------------------------------------------------------------------
# read in known metabolites shared by 2 datasets
metMap <- read.table(metMapUniqFile,sep="\t",header=T)
metMap <- metMap[complete.cases(metMap),]
colnames(metMap) <- c("Metabolite1","Metabolite2")
dim(metMap)

# read in Dataset1
metData1 <- read.table(metDataFile1,sep="\t",header=T,row.names=1,check.names=F)
dim(metData1)

knownData1 <- metData1[,match(metMap$Metabolite1,colnames(metData1))]
dim(knownData1)

# read in Dataset2
metData2 <- read.table(metDataFile2,sep="\t",header=T,row.names=1,check.names=F)
dim(metData2)

knownData2 <- metData2[,match(metMap$Metabolite2,colnames(metData2))]
dim(knownData2)

# use column headers for knownData1 (for doing lm below)
colnames(knownData2) <- colnames(knownData1)

# ----------------------------------------------------------------------------
# IMPUTATION

# impute Dataset2 signal values for Dataset1 samples
impData1 <- NULL
for (m in 1:dim(metData2)[2]) {
	if (colnames(metData2)[m] %in% metMap$Metabolite2) { # for shared known metabolites
		n <- which(metMap$Metabolite2 == colnames(metData2)[m])
		impData1 <- cbind(impData1, predict(lm(metData2[,m] ~ ., data=knownData2[,-n]), knownData1[,-n]))
	
	} else { # for all other signals
		impData1 <- cbind(impData1, predict(lm(metData2[,m] ~ ., data=knownData2), knownData1))
		
	}
}

colnames(impData1) <- colnames(metData2)

# output imputed data
write.table(impData1,out_impFile1,sep="\t",quote=F)


# impute Dataset1 signal values for Dataset2 samples
impData2 <- NULL
for (m in 1:dim(metData1)[2]) {
	if (colnames(metData1)[m] %in% metMap$Metabolite1) { # for shared known metabolites
		n <- which(metMap$Metabolite1 == colnames(metData1)[m])
		impData2 <- cbind(impData2, predict(lm(metData1[,m] ~ ., data=knownData1[,-n]), knownData2[,-n]))
	
	} else { # for all other signals
		impData2 <- cbind(impData2, predict(lm(metData1[,m] ~ ., data=knownData1), knownData2))
		
	}
}

colnames(impData2) <- colnames(metData1)

# output imputed data
write.table(impData2,out_impFile2,sep="\t",quote=F)

# PRINT PROCESS TIME
print("TOTAL SCRIPT TIME:")
print(proc.time() - ptm)

