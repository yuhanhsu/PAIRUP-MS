##########################################################################################
### Match signals across datasets using m/z + imputed correlation
### generate 6 sets of matches:
###	matching Dataset1 signals to Dataset2 or vice versa
### 	calculate correlation using Dataset1, Dataset2, or all samples
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# TIME PROCESS
ptm <- proc.time()

# read in parameters from input argument file name
source(commandArgs(trailingOnly=T)[1],echo=T)


# ----------------------------------------------------------------------------------
# read in observed datasets
print("OBSERVED DATA:")

metData1 <- read.table(metDataFile1,sep="\t",header=T,row.names=1,check.names=F)
dim(metData1)

metData2 <- read.table(metDataFile2,sep="\t",header=T,row.names=1,check.names=F)
dim(metData2)

# read in imputed datasets
print("IMPUTED DATA:")

impData1 <- read.table(out_impFile1,sep="\t",header=T,row.names=1,check.names=F)
dim(impData1)

impData2 <- read.table(out_impFile2,sep="\t",header=T,row.names=1,check.names=F)
dim(impData2)

# read in m/z and RT information
print("METABOLITES WITH M/Z and RT:")

metInfo1 <- NULL
for (i in 1:length(methods1)) {
	tempTable <- read.table(metInfoFiles1[i],sep="\t",header=T,stringsAsFactors=F)
	tempTable$Metabolite <- gsub(" ","_",tempTable$Metabolite)
	tempTable$Method <- methods1[i]
	tempTable$IonMode <- ionModes1[i]
	metInfo1 <- rbind(metInfo1,tempTable)
}
metInfo1 <- metInfo1[complete.cases(metInfo1),]
metInfo1 <- metInfo1[metInfo1$Metabolite %in% colnames(metData1),]
print(dim(metInfo1))

metInfo2 <- NULL
for (i in 1:length(methods2)) {
	tempTable <- read.table(metInfoFiles2[i],sep="\t",header=T,stringsAsFactors=F)
	
	tempName <- gsub(" ","_",tempTable$Metabolite)
	tempName2 <- paste(tempName,"_",methods2[i],sep="")

	tempTable$Metabolite <- NA
	tempTable$Metabolite[tempName %in% colnames(metData2)] <- tempName[tempName %in% colnames(metData2)]
	tempTable$Metabolite[tempName2 %in% colnames(metData2)] <- tempName2[tempName2 %in% colnames(metData2)]
	
	tempTable$Method <- methods2[i]
	tempTable$IonMode <- ionModes2[i]
	metInfo2 <- rbind(metInfo2,tempTable)
}
metInfo2 <- metInfo2[complete.cases(metInfo2),]
metInfo2 <- metInfo2[metInfo2$Metabolite %in% colnames(metData2),]
print(dim(metInfo2))


# ---------------------------------------------------------------------------
# MATCHING: match signals using m/z + correlation between observed and imputed data
source("../src/Match_MZ_Function.r",echo=T)
vMatchMZ <- Vectorize(matchMZ) # vectorize matchMZ function

print("MATCH 1-2 USING IMP 1...")
match12_imp1_Table <- NULL 

if (dim(metInfo1)[1] > 0 & dim(metInfo2)[1] > 0) {
	for (i in 1:dim(metInfo1)[1]) {
		isMzMatch <- vMatchMZ(metInfo1$m.z[i],metInfo2$m.z,metInfo1$IonMode[i],metInfo2$IonMode,mzMode,mzWindow)

		if (sum(isMzMatch) > 0) {
			mets <- metInfo2$Metabolite[isMzMatch]
			corrs <- cor(metData1[,metInfo1$Metabolite[i]],impData1[,mets])
			rankedCorrs <- corrs[order(-corrs)]
			rankedMets <- mets[order(-corrs)]

			m <- which(metInfo2$Metabolite==rankedMets[1]) # match!

			tempRow <- c(metInfo1[i,c("Metabolite","Method","IonMode","m.z","RT")],
					metInfo2[m,c("Metabolite","Method","IonMode","m.z","RT")],
					DeltaMZ=metInfo1[i,"m.z"]-metInfo2[m,"m.z"],DeltaRT=metInfo1[i,"RT"]-metInfo2[m,"RT"],
					Match_Corr=rankedCorrs[1])
			match12_imp1_Table <- rbind(match12_imp1_Table,tempRow)
		}
		if (i %% 1000 == 0) {
			print(i)
		}
	}
}

print("MATCH 1-2 USING IMP 1: # OF MATCHED PAIRS:")
print(dim(match12_imp1_Table)[1])

# output matched pairs
colnames(match12_imp1_Table) <- c("Metabolite1","Method1","Ion1","MZ1","RT1","Metabolite2","Method2","Ion2","MZ2","RT2",
	"DeltaMZ","DeltaRT","Match_Corr")
write.table(match12_imp1_Table,file=out_match12_imp1_File,sep="\t",quote=F,row.names=F)


# ---------------------------------------
print("MATCH 2-1 USING IMP 1...")
match21_imp1_Table <- NULL 

if (dim(metInfo1)[1] > 0 & dim(metInfo2)[1] > 0) {
	for (i in 1:dim(metInfo2)[1]) {
		isMzMatch <- vMatchMZ(metInfo2$m.z[i],metInfo1$m.z,metInfo2$IonMode[i],metInfo1$IonMode,mzMode,mzWindow)

		if (sum(isMzMatch) > 0) {
			mets <- metInfo1$Metabolite[isMzMatch]
			corrs <- cor(impData1[,metInfo2$Metabolite[i]],metData1[,mets])
			rankedCorrs <- corrs[order(-corrs)]
			rankedMets <- mets[order(-corrs)]

			m <- which(metInfo1$Metabolite==rankedMets[1]) # match!

			tempRow <- c(metInfo2[i,c("Metabolite","Method","IonMode","m.z","RT")],
					metInfo1[m,c("Metabolite","Method","IonMode","m.z","RT")],
					DeltaMZ=metInfo2[i,"m.z"]-metInfo1[m,"m.z"],DeltaRT=metInfo2[i,"RT"]-metInfo1[m,"RT"],
					Match_Corr=rankedCorrs[1])
			match21_imp1_Table <- rbind(match21_imp1_Table,tempRow)
		}
		if (i %% 1000 == 0) {
			print(i)
		}
	}
}

print("MATCH 2-1 USING IMP 1: # OF MATCHED PAIRS:")
print(dim(match21_imp1_Table)[1])

# output matched pairs
colnames(match21_imp1_Table) <- c("Metabolite1","Method1","Ion1","MZ1","RT1","Metabolite2","Method2","Ion2","MZ2","RT2",
	"DeltaMZ","DeltaRT","Match_Corr")
write.table(match21_imp1_Table,file=out_match21_imp1_File,sep="\t",quote=F,row.names=F)


# ---------------------------------------
print("MATCH 1-2 USING IMP 2...")
match12_imp2_Table <- NULL 

if (dim(metInfo1)[1] > 0 & dim(metInfo2)[1] > 0) {
	for (i in 1:dim(metInfo1)[1]) {
		isMzMatch <- vMatchMZ(metInfo1$m.z[i],metInfo2$m.z,metInfo1$IonMode[i],metInfo2$IonMode,mzMode,mzWindow)

		if (sum(isMzMatch) > 0) {
			mets <- metInfo2$Metabolite[isMzMatch]
			corrs <- cor(impData2[,metInfo1$Metabolite[i]],metData2[,mets])
			rankedCorrs <- corrs[order(-corrs)]
			rankedMets <- mets[order(-corrs)]

			m <- which(metInfo2$Metabolite==rankedMets[1]) # match!

			tempRow <- c(metInfo1[i,c("Metabolite","Method","IonMode","m.z","RT")],
					metInfo2[m,c("Metabolite","Method","IonMode","m.z","RT")],
					DeltaMZ=metInfo1[i,"m.z"]-metInfo2[m,"m.z"],DeltaRT=metInfo1[i,"RT"]-metInfo2[m,"RT"],
					Match_Corr=rankedCorrs[1])
			match12_imp2_Table <- rbind(match12_imp2_Table,tempRow)
		}
		if (i %% 1000 == 0) {
			print(i)
		}
	}
}

print("MATCH 1-2 USING IMP 2: # OF MATCHED PAIRS:")
print(dim(match12_imp2_Table)[1])

# output matched pairs
colnames(match12_imp2_Table) <- c("Metabolite1","Method1","Ion1","MZ1","RT1","Metabolite2","Method2","Ion2","MZ2","RT2",
	"DeltaMZ","DeltaRT","Match_Corr")
write.table(match12_imp2_Table,file=out_match12_imp2_File,sep="\t",quote=F,row.names=F)


# ---------------------------------------
print("MATCH 2-1 USING IMP 2...")
match21_imp2_Table <- NULL 

if (dim(metInfo1)[1] > 0 & dim(metInfo2)[1] > 0) {
	for (i in 1:dim(metInfo2)[1]) {
		isMzMatch <- vMatchMZ(metInfo2$m.z[i],metInfo1$m.z,metInfo2$IonMode[i],metInfo1$IonMode,mzMode,mzWindow)

		if (sum(isMzMatch) > 0) {
			mets <- metInfo1$Metabolite[isMzMatch]
			corrs <- cor(metData2[,metInfo2$Metabolite[i]],impData2[,mets])
			rankedCorrs <- corrs[order(-corrs)]
			rankedMets <- mets[order(-corrs)]

			m <- which(metInfo1$Metabolite==rankedMets[1]) # match!

			tempRow <- c(metInfo2[i,c("Metabolite","Method","IonMode","m.z","RT")],
					metInfo1[m,c("Metabolite","Method","IonMode","m.z","RT")],
					DeltaMZ=metInfo2[i,"m.z"]-metInfo1[m,"m.z"],DeltaRT=metInfo2[i,"RT"]-metInfo1[m,"RT"],
					Match_Corr=rankedCorrs[1])
			match21_imp2_Table <- rbind(match21_imp2_Table,tempRow)
		}
		if (i %% 1000 == 0) {
			print(i)
		}
	}
}

print("MATCH 2-1 USING IMP 2: # OF MATCHED PAIRS:")
print(dim(match21_imp2_Table)[1])

# output matched pairs
colnames(match21_imp2_Table) <- c("Metabolite1","Method1","Ion1","MZ1","RT1","Metabolite2","Method2","Ion2","MZ2","RT2",
	"DeltaMZ","DeltaRT","Match_Corr")
write.table(match21_imp2_Table,file=out_match21_imp2_File,sep="\t",quote=F,row.names=F)


# ---------------------------------------
print("MATCH 1-2 USING IMP 1+2...") # correlation across all samples in both datasets
match12_impAll_Table <- NULL 

if (dim(metInfo1)[1] > 0 & dim(metInfo2)[1] > 0) {
	for (i in 1:dim(metInfo1)[1]) {
		isMzMatch <- vMatchMZ(metInfo1$m.z[i],metInfo2$m.z,metInfo1$IonMode[i],metInfo2$IonMode,mzMode,mzWindow)

		if (sum(isMzMatch) > 0) {
			mets <- metInfo2$Metabolite[isMzMatch]

			if (sum(isMzMatch) > 1) {
				corrs <- cor(c(metData1[,metInfo1$Metabolite[i]],impData2[,metInfo1$Metabolite[i]]),
						rbind(impData1[,mets],metData2[,mets]))
			} else {
				corrs <- cor(c(metData1[,metInfo1$Metabolite[i]],impData2[,metInfo1$Metabolite[i]]),
						c(impData1[,mets],metData2[,mets]))
			}
			rankedCorrs <- corrs[order(-corrs)]
			rankedMets <- mets[order(-corrs)]

			m <- which(metInfo2$Metabolite==rankedMets[1]) # match!

			tempRow <- c(metInfo1[i,c("Metabolite","Method","IonMode","m.z","RT")],
					metInfo2[m,c("Metabolite","Method","IonMode","m.z","RT")],
					DeltaMZ=metInfo1[i,"m.z"]-metInfo2[m,"m.z"],DeltaRT=metInfo1[i,"RT"]-metInfo2[m,"RT"],
					Match_Corr=rankedCorrs[1])
			match12_impAll_Table <- rbind(match12_impAll_Table,tempRow)
		}
		if (i %% 1000 == 0) {
			print(i)
		}
	}
}

print("MATCH 1-2 USING IMP 1+2: # OF MATCHED PAIRS:")
print(dim(match12_impAll_Table)[1])

# output matched pairs
colnames(match12_impAll_Table) <- c("Metabolite1","Method1","Ion1","MZ1","RT1","Metabolite2","Method2","Ion2","MZ2","RT2",
	"DeltaMZ","DeltaRT","Match_Corr")
write.table(match12_impAll_Table,file=out_match12_impAll_File,sep="\t",quote=F,row.names=F)


# ---------------------------------------
print("MATCH 2-1 USING IMP 1+2...")
match21_impAll_Table <- NULL 

if (dim(metInfo1)[1] > 0 & dim(metInfo2)[1] > 0) {
	for (i in 1:dim(metInfo2)[1]) {
		isMzMatch <- vMatchMZ(metInfo2$m.z[i],metInfo1$m.z,metInfo2$IonMode[i],metInfo1$IonMode,mzMode,mzWindow)

		if (sum(isMzMatch) > 0) {
			mets <- metInfo1$Metabolite[isMzMatch]
 
			if (sum(isMzMatch) > 1) {
				corrs <- cor(c(metData2[,metInfo2$Metabolite[i]],impData1[,metInfo2$Metabolite[i]]),
						rbind(impData2[,mets],metData1[,mets]))
			} else {
				corrs <- cor(c(metData2[,metInfo2$Metabolite[i]],impData1[,metInfo2$Metabolite[i]]),
						c(impData2[,mets],metData1[,mets]))
			}
				
			rankedCorrs <- corrs[order(-corrs)]
			rankedMets <- mets[order(-corrs)]

			m <- which(metInfo1$Metabolite==rankedMets[1]) # match!

			tempRow <- c(metInfo2[i,c("Metabolite","Method","IonMode","m.z","RT")],
					metInfo1[m,c("Metabolite","Method","IonMode","m.z","RT")],
					DeltaMZ=metInfo2[i,"m.z"]-metInfo1[m,"m.z"],DeltaRT=metInfo2[i,"RT"]-metInfo1[m,"RT"],
					Match_Corr=rankedCorrs[1])
			match21_impAll_Table <- rbind(match21_impAll_Table,tempRow)
		}
		if (i %% 1000 == 0) {
			print(i)
		}
	}
}

print("MATCH 2-1 USING IMP 1+2: # OF MATCHED PAIRS:")
print(dim(match21_impAll_Table)[1])

# output matched pairs
colnames(match21_impAll_Table) <- c("Metabolite1","Method1","Ion1","MZ1","RT1","Metabolite2","Method2","Ion2","MZ2","RT2",
	"DeltaMZ","DeltaRT","Match_Corr")
write.table(match21_impAll_Table,file=out_match21_impAll_File,sep="\t",quote=F,row.names=F)


# PRINT PROCESS TIME
print("TOTAL SCRIPT TIME:")
print(proc.time() - ptm)

