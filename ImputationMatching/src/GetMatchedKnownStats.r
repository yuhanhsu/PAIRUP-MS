##########################################################################################
### Get calibration statistics using matched shared known metabolites
###
### author: Yu-Han Hsu 
##########################################################################################

rm(list=ls())

library(ggplot2)

# read in parameters from input argument file name
source(commandArgs(trailingOnly=T)[1],echo=T)


# ----------------------------------------------------------------------------------
# read in known metabolites shared by 2 datasets
metMap <- read.table(metMapFile,sep="\t",header=T,stringsAsFactors=F)
metMap <- metMap[complete.cases(metMap),]
colnames(metMap) <- c("Metabolite1","Metabolite2")
dim(metMap)

# read in observed datasets
metData1 <- read.table(metDataFile1,sep="\t",header=T,row.names=1,check.names=F)
dim(metData1)

metData2 <- read.table(metDataFile2,sep="\t",header=T,row.names=1,check.names=F)
dim(metData2)

# read in imputed datasets
impData1 <- read.table(out_impFile1,sep="\t",header=T,row.names=1,check.names=F)
dim(impData1)

impData2 <- read.table(out_impFile2,sep="\t",header=T,row.names=1,check.names=F)
dim(impData2)

# read in matched pairs
match12_rtTable <- read.table(out_match12_RtFile,sep="\t",header=T,stringsAsFactors=F)
match21_rtTable <- read.table(out_match21_RtFile,sep="\t",header=T,stringsAsFactors=F)

match12_imp1_Table <- read.table(out_match12_imp1_File,sep="\t",header=T,stringsAsFactors=F)
match21_imp1_Table <- read.table(out_match21_imp1_File,sep="\t",header=T,stringsAsFactors=F)
match12_imp2_Table <- read.table(out_match12_imp2_File,sep="\t",header=T,stringsAsFactors=F)
match21_imp2_Table <- read.table(out_match21_imp2_File,sep="\t",header=T,stringsAsFactors=F)
match12_impAll_Table <- read.table(out_match12_impAll_File,sep="\t",header=T,stringsAsFactors=F)
match21_impAll_Table <- read.table(out_match21_impAll_File,sep="\t",header=T,stringsAsFactors=F)

# -----------------------------------------------------------------------------------------------
# Get Shared Knowns Matching Results

# MATCH 1-2 USING RT
knownMatch12_rtTable <- match12_rtTable[match12_rtTable$Metabolite1 %in% metMap$Metabolite1,]
#dim(knownMatch12_rtTable)[1]

knownMatch12_rtTable$RealMetabolite <- NA
knownMatch12_rtTable$CorrectMatch <- FALSE
knownMatch12_rtTable$Match_Real_ObsCorr <- 0 # correlation between match and real in observed data
for (i in 1:dim(knownMatch12_rtTable)[1]) {
	met <- knownMatch12_rtTable$Metabolite1[i]
	match <- knownMatch12_rtTable$Metabolite2[i]
	real <- metMap$Metabolite2[metMap$Metabolite1==met]

	knownMatch12_rtTable$RealMetabolite[i] <- paste(real,collapse=";")

	if (match %in% real) {
		knownMatch12_rtTable$CorrectMatch[i] <- TRUE
	}

	knownMatch12_rtTable$Match_Real_ObsCorr[i] <- max(cor(metData2[,match],metData2[,real]))
}
write.table(format(knownMatch12_rtTable,digits=4),out_knownMatch12_rtFile,sep="\t",quote=F,row.names=F)


# MATCH 2-1 USING RT
knownMatch21_rtTable <- match21_rtTable[match21_rtTable$Metabolite1 %in% metMap$Metabolite2,]
#dim(knownMatch21_rtTable)[1]

knownMatch21_rtTable$RealMetabolite <- NA
knownMatch21_rtTable$CorrectMatch <- FALSE
knownMatch21_rtTable$Match_Real_ObsCorr <- 0 # correlation between match and real in observed data
for (i in 1:dim(knownMatch21_rtTable)[1]) {
	met <- knownMatch21_rtTable$Metabolite1[i]
	match <- knownMatch21_rtTable$Metabolite2[i]
	real <- metMap$Metabolite1[metMap$Metabolite2==met]

	knownMatch21_rtTable$RealMetabolite[i] <- paste(real,collapse=";")

	if (match %in% real) {
		knownMatch21_rtTable$CorrectMatch[i] <- TRUE
	}

	knownMatch21_rtTable$Match_Real_ObsCorr[i] <- max(cor(metData1[,match],metData1[,real]))
}
write.table(format(knownMatch21_rtTable,digits=4),out_knownMatch21_rtFile,sep="\t",quote=F,row.names=F)


# MATCH 1-2 USING IMP 1
knownMatch12_imp1_Table <- match12_imp1_Table[match12_imp1_Table$Metabolite1 %in% metMap$Metabolite1,]
#dim(knownMatch12_imp1_Table)[1]

knownMatch12_imp1_Table$RealMetabolite <- NA
knownMatch12_imp1_Table$CorrectMatch <- FALSE
knownMatch12_imp1_Table$Real_Corr <- 0 # correlation between known and real in imputed data
knownMatch12_imp1_Table$Match_Real_ObsCorr <- 0 # correlation between match and real in observed data
for (i in 1:dim(knownMatch12_imp1_Table)[1]) {
	met <- knownMatch12_imp1_Table$Metabolite1[i]
	match <- knownMatch12_imp1_Table$Metabolite2[i]
	real <- metMap$Metabolite2[metMap$Metabolite1==met]

	knownMatch12_imp1_Table$RealMetabolite[i] <- paste(real,collapse=";")

	if (match %in% real) {
		knownMatch12_imp1_Table$CorrectMatch[i] <- TRUE
	}

	knownMatch12_imp1_Table$Real_Corr[i] <- max(cor(metData1[,met],impData1[,real]))
	knownMatch12_imp1_Table$Match_Real_ObsCorr[i] <- max(cor(metData2[,match],metData2[,real]))
}
write.table(format(knownMatch12_imp1_Table,digits=4),out_knownMatch12_imp1_File,sep="\t",quote=F,row.names=F)


# MATCH 2-1 USING IMP 1
knownMatch21_imp1_Table <- match21_imp1_Table[match21_imp1_Table$Metabolite1 %in% metMap$Metabolite2,]
#dim(knownMatch21_imp1_Table)[1]

knownMatch21_imp1_Table$RealMetabolite <- NA
knownMatch21_imp1_Table$CorrectMatch <- FALSE
knownMatch21_imp1_Table$Real_Corr <- 0 # correlation between known and real in imputed data
knownMatch21_imp1_Table$Match_Real_ObsCorr <- 0 # correlation between match and real in observed data
for (i in 1:dim(knownMatch21_imp1_Table)[1]) {
	met <- knownMatch21_imp1_Table$Metabolite1[i]
	match <- knownMatch21_imp1_Table$Metabolite2[i]
	real <- metMap$Metabolite1[metMap$Metabolite2==met]

	knownMatch21_imp1_Table$RealMetabolite[i] <- paste(real,collapse=";")

	if (match %in% real) {
		knownMatch21_imp1_Table$CorrectMatch[i] <- TRUE
	}

	knownMatch21_imp1_Table$Real_Corr[i] <- max(cor(impData1[,met],metData1[,real]))
	knownMatch21_imp1_Table$Match_Real_ObsCorr[i] <- max(cor(metData1[,match],metData1[,real]))
}
write.table(format(knownMatch21_imp1_Table,digits=4),out_knownMatch21_imp1_File,sep="\t",quote=F,row.names=F)


# MATCH 1-2 USING IMP 2
knownMatch12_imp2_Table <- match12_imp2_Table[match12_imp2_Table$Metabolite1 %in% metMap$Metabolite1,]
#dim(knownMatch12_imp2_Table)[1]

knownMatch12_imp2_Table$RealMetabolite <- NA
knownMatch12_imp2_Table$CorrectMatch <- FALSE
knownMatch12_imp2_Table$Real_Corr <- 0 # correlation between known and real in imputed data
knownMatch12_imp2_Table$Match_Real_ObsCorr <- 0 # correlation between match and real in observed data
for (i in 1:dim(knownMatch12_imp2_Table)[1]) {
	met <- knownMatch12_imp2_Table$Metabolite1[i]
	match <- knownMatch12_imp2_Table$Metabolite2[i]
	real <- metMap$Metabolite2[metMap$Metabolite1==met]

	knownMatch12_imp2_Table$RealMetabolite[i] <- paste(real,collapse=";")

	if (match %in% real) {
		knownMatch12_imp2_Table$CorrectMatch[i] <- TRUE
	}

	knownMatch12_imp2_Table$Real_Corr[i] <- max(cor(impData2[,met],metData2[,real]))
	knownMatch12_imp2_Table$Match_Real_ObsCorr[i] <- max(cor(metData2[,match],metData2[,real]))
}
write.table(format(knownMatch12_imp2_Table,digits=4),out_knownMatch12_imp2_File,sep="\t",quote=F,row.names=F)


# MATCH 2-1 USING IMP 2
knownMatch21_imp2_Table <- match21_imp2_Table[match21_imp2_Table$Metabolite1 %in% metMap$Metabolite2,]
#dim(knownMatch21_imp2_Table)[1]

knownMatch21_imp2_Table$RealMetabolite <- NA
knownMatch21_imp2_Table$CorrectMatch <- FALSE
knownMatch21_imp2_Table$Real_Corr <- 0 # correlation between known and real in imputed data
knownMatch21_imp2_Table$Match_Real_ObsCorr <- 0 # correlation between match and real in observed data
for (i in 1:dim(knownMatch21_imp2_Table)[1]) {
	met <- knownMatch21_imp2_Table$Metabolite1[i]
	match <- knownMatch21_imp2_Table$Metabolite2[i]
	real <- metMap$Metabolite1[metMap$Metabolite2==met]

	knownMatch21_imp2_Table$RealMetabolite[i] <- paste(real,collapse=";")

	if (match %in% real) {
		knownMatch21_imp2_Table$CorrectMatch[i] <- TRUE
	}

	knownMatch21_imp2_Table$Real_Corr[i] <- max(cor(metData2[,met],impData2[,real]))
	knownMatch21_imp2_Table$Match_Real_ObsCorr[i] <- max(cor(metData1[,match],metData1[,real]))
}
write.table(format(knownMatch21_imp2_Table,digits=4),out_knownMatch21_imp2_File,sep="\t",quote=F,row.names=F)


# MATCH 1-2 USING IMP 1+2
knownMatch12_impAll_Table <- match12_impAll_Table[match12_impAll_Table$Metabolite1 %in% metMap$Metabolite1,]
#dim(knownMatch12_impAll_Table)[1]

knownMatch12_impAll_Table$RealMetabolite <- NA
knownMatch12_impAll_Table$CorrectMatch <- FALSE
knownMatch12_impAll_Table$Real_Corr <- 0 # correlation between known and real in imputed data
knownMatch12_impAll_Table$Match_Real_ObsCorr <- 0 # correlation between match and real in observed data
for (i in 1:dim(knownMatch12_impAll_Table)[1]) {
	met <- knownMatch12_impAll_Table$Metabolite1[i]
	match <- knownMatch12_impAll_Table$Metabolite2[i]
	real <- metMap$Metabolite2[metMap$Metabolite1==met]

	knownMatch12_impAll_Table$RealMetabolite[i] <- paste(real,collapse=";")

	if (match %in% real) {
		knownMatch12_impAll_Table$CorrectMatch[i] <- TRUE
	}

	if (length(real) > 1) {
		knownMatch12_impAll_Table$Real_Corr[i] <- max(cor(c(metData1[,met],impData2[,met]),
								rbind(impData1[,real],metData2[,real])))
	} else {
		knownMatch12_impAll_Table$Real_Corr[i] <- cor(c(metData1[,met],impData2[,met]),
								c(impData1[,real],metData2[,real]))
	}
	
	knownMatch12_impAll_Table$Match_Real_ObsCorr[i] <- max(cor(metData2[,match],metData2[,real]))
}
write.table(format(knownMatch12_impAll_Table,digits=4),out_knownMatch12_impAll_File,sep="\t",quote=F,row.names=F)


# MATCH 2-1 USING IMP 1+2
knownMatch21_impAll_Table <- match21_impAll_Table[match21_impAll_Table$Metabolite1 %in% metMap$Metabolite2,]
#dim(knownMatch21_impAll_Table)[1]

knownMatch21_impAll_Table$RealMetabolite <- NA
knownMatch21_impAll_Table$CorrectMatch <- FALSE
knownMatch21_impAll_Table$Real_Corr <- 0 # correlation between known and real in imputed data
knownMatch21_impAll_Table$Match_Real_ObsCorr <- 0 # correlation between match and real in observed data
for (i in 1:dim(knownMatch21_impAll_Table)[1]) {
	met <- knownMatch21_impAll_Table$Metabolite1[i]
	match <- knownMatch21_impAll_Table$Metabolite2[i]
	real <- metMap$Metabolite1[metMap$Metabolite2==met]

	knownMatch21_impAll_Table$RealMetabolite[i] <- paste(real,collapse=";")

	if (match %in% real) {
		knownMatch21_impAll_Table$CorrectMatch[i] <- TRUE
	}

	if (length(real) > 1) {
		knownMatch21_impAll_Table$Real_Corr[i] <- max(cor(c(impData1[,met],metData2[,met]),
								rbind(metData1[,real],impData2[,real])))
	} else {
		knownMatch21_impAll_Table$Real_Corr[i] <- cor(c(impData1[,met],metData2[,met]),
								c(metData1[,real],impData2[,real]))
	}
	
	knownMatch21_impAll_Table$Match_Real_ObsCorr[i] <- max(cor(metData1[,match],metData1[,real]))
}
write.table(format(knownMatch21_impAll_Table,digits=4),out_knownMatch21_impAll_File,sep="\t",quote=F,row.names=F)


# -----------------------------------------------------------------------------------------------
# PLOTS: only for IMP-matched tables

tableList <- list(knownMatch12_imp1_Table,knownMatch21_imp1_Table,knownMatch12_imp2_Table,knownMatch21_imp2_Table,
			knownMatch12_impAll_Table,knownMatch21_impAll_Table)
plotList <- list(out_knownMatch12_imp1_plotFile,out_knownMatch21_imp1_plotFile,
			out_knownMatch12_imp2_plotFile,out_knownMatch21_imp2_plotFile,
			out_knownMatch12_impAll_plotFile,out_knownMatch21_impAll_plotFile)

# compare correlation between metabolite, match, and real
for (i in 1:length(tableList)) {
	tempTable <- tableList[[i]]
	tempPlot <- plotList[[i]]
	
	pdf(tempPlot,width=6,height=4)
	hist(tempTable$Match_Corr ^ 2,main="Known-Match R^2 in Imputed Data")
	hist(tempTable$Real_Corr ^ 2,main="Known-Real R^2 in Imputed Data")
	hist(tempTable$Match_Real_ObsCorr ^ 2,main="Match-Real R^2 in Observed Data")

	print(ggplot(tempTable,aes(x=Match_Corr^2,y=Real_Corr^2,color=CorrectMatch)) + geom_point() + 
		xlab("Known-Match R^2 in Imputed Data") + ylab("Known-Real R^2 in Imputed Data") + theme_bw())

	print(ggplot(tempTable,aes(x=Match_Corr^2,y=Match_Real_ObsCorr^2,color=CorrectMatch)) + geom_point() + 
		xlab("Known-Match R^2 in Imputed Data") + ylab("Match-Real R^2 in Observed Data") + theme_bw())

	dev.off()
}

