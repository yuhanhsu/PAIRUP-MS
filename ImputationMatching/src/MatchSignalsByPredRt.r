##########################################################################################
### Match signals across datasets using m/z + predicted RT
### (use shared knowns to build linear model for predicting RT)
### generate 2 sets of matches:
###     matching Dataset1 signals to Dataset2 or vice versa
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

# TIME PROCESS
ptm <- proc.time()

# read in parameters from input argument file name
source(commandArgs(trailingOnly=T)[1],echo=T)


# ----------------------------------------------------------------------------------
# read in metabolite signals included in each dataset
print("METABOLITE LISTS:")

metList1 <- read.table(metListFile1,header=F,stringsAsFactors=F)$V1
print(length(metList1))

metList2 <- read.table(metListFile2,header=F,stringsAsFactors=F)$V1
print(length(metList2))

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
metInfo1 <- metInfo1[metInfo1$Metabolite %in% metList1,]
print(dim(metInfo1))

metInfo2 <- NULL
for (i in 1:length(methods2)) {
	tempTable <- read.table(metInfoFiles2[i],sep="\t",header=T,stringsAsFactors=F)

	tempName <- gsub(" ","_",tempTable$Metabolite)
	tempName2 <- paste(tempName,"_",methods2[i],sep="")
	
	tempTable$Metabolite <- NA
	tempTable$Metabolite[tempName %in% metList2] <- tempName[tempName %in% metList2]
	tempTable$Metabolite[tempName2 %in% metList2] <- tempName2[tempName2 %in% metList2]

	tempTable$Method <- methods2[i]
	tempTable$IonMode <- ionModes2[i]
	metInfo2 <- rbind(metInfo2,tempTable)
}
metInfo2 <- metInfo2[complete.cases(metInfo2),]
metInfo2 <- metInfo2[metInfo2$Metabolite %in% metList2,]
print(dim(metInfo2))

# ---------------------------------------------------------------------------
# BUILD RT PREDICTON MODELS USING SHARED KNOWNS
knownMapTable <- read.table(metMapFile,header=T)
knownMapTable <- knownMapTable[complete.cases(knownMapTable),]
print(dim(knownMapTable))

knownInfo1 <- metInfo1[match(knownMapTable[,1],metInfo1$Metabolite),c("Metabolite","Method","m.z","RT")]
knownInfo2 <- metInfo2[match(knownMapTable[,2],metInfo2$Metabolite),c("Metabolite","Method","m.z","RT")]
sharedKnownInfo <- cbind(knownInfo1,knownInfo2)
sharedKnownInfo <- sharedKnownInfo[complete.cases(sharedKnownInfo),]
colnames(sharedKnownInfo) <- c("Metabolite1","Method1","m.z1","RT1","Metabolite2","Method2","m.z2","RT2")
print(dim(sharedKnownInfo))

model12 <- lm(RT2 ~ RT1,data=sharedKnownInfo)
print(summary(model12))

model21 <- lm(RT1 ~ RT2,data=sharedKnownInfo)
print(summary(model21))

# ---------------------------------------------------------------------------
# MATCHING: match signals using m/z + closeness to predicted RT 
source("../src/Match_MZ_Function.r",echo=T)
vMatchMZ <- Vectorize(matchMZ) # vectorize matchMZ function

print("MATCH 1-2")
match12_Table <- NULL 

if (dim(metInfo1)[1] > 0 & dim(metInfo2)[1] > 0) {
	for (i in 1:dim(metInfo1)[1]) {
		isMzMatch <- vMatchMZ(metInfo1$m.z[i],metInfo2$m.z,metInfo1$IonMode[i],metInfo2$IonMode,mzMode,mzWindow)

		if (sum(isMzMatch) > 0) {
			# calculate predicted RT using linear model
			predRT <- predict(model12,data.frame(RT1=metInfo1$RT[i]))
			
			# calculate distance from predicted RT
			deltaRT <- abs(metInfo2$RT[isMzMatch] - predRT)
			rankedMet2 <- metInfo2$Metabolite[isMzMatch][order(deltaRT)] 

			m <- which(metInfo2$Metabolite==rankedMet2[1]) # match!

			tempRow <- data.frame(metInfo1[i,c("Metabolite","Method","IonMode","m.z","RT")],
					metInfo2[m,c("Metabolite","Method","IonMode","m.z","RT")],
					DeltaMZ=round(metInfo1[i,"m.z"]-metInfo2[m,"m.z"],4),
					PredRT=round(predRT,4),DeltaPredRT=round(abs(metInfo2[m,"RT"]-predRT),4))
			match12_Table <- rbind(match12_Table,tempRow)
		}
		if (i %% 1000 == 0) {
			print(i)
		}
	}
}

print("# OF MATCHED PAIRS:")
print(dim(match12_Table)[1])

# output matched pairs
colnames(match12_Table) <- c("Metabolite1","Method1","Ion1","MZ1","RT1","Metabolite2","Method2","Ion2","MZ2","RT2",
	"DeltaMZ","PredRT","DeltaPredRT")
write.table(match12_Table,file=out_match12_RtFile,sep="\t",quote=F,row.names=F)

# ----------------
print("MATCH 2-1")
match21_Table <- NULL 

if (dim(metInfo1)[1] > 0 & dim(metInfo2)[1] > 0) {
	for (i in 1:dim(metInfo2)[1]) {
		isMzMatch <- vMatchMZ(metInfo2$m.z[i],metInfo1$m.z,metInfo2$IonMode[i],metInfo1$IonMode,mzMode,mzWindow)

		if (sum(isMzMatch) > 0) {
			# calculate predicted RT using linear model
			predRT <- predict(model21,data.frame(RT2=metInfo2$RT[i]))
			
			# calculate distance from predicted RT
			deltaRT <- abs(metInfo1$RT[isMzMatch] - predRT)
			rankedMet1 <- metInfo1$Metabolite[isMzMatch][order(deltaRT)] 

			m <- which(metInfo1$Metabolite==rankedMet1[1]) # match!

			tempRow <- data.frame(metInfo2[i,c("Metabolite","Method","IonMode","m.z","RT")],
					metInfo1[m,c("Metabolite","Method","IonMode","m.z","RT")],
					DeltaMZ=round(metInfo2[i,"m.z"]-metInfo1[m,"m.z"],4),
                                        PredRT=round(predRT,4),DeltaPredRT=round(abs(metInfo1[m,"RT"]-predRT),4))
			match21_Table <- rbind(match21_Table,tempRow)
		}
		if (i %% 1000 == 0) {
			print(i)
		}
	}
}

print("# OF MATCHED PAIRS:")
print(dim(match21_Table)[1])

# output matched pairs
colnames(match21_Table) <- c("Metabolite1","Method1","Ion1","MZ1","RT1","Metabolite2","Method2","Ion2","MZ2","RT2",
	"DeltaMZ","PredRT","DeltaPredRT")
write.table(match21_Table,file=out_match21_RtFile,sep="\t",quote=F,row.names=F)


# PRINT PROCESS TIME
print("TOTAL SCRIPT TIME:")
print(proc.time() - ptm)

