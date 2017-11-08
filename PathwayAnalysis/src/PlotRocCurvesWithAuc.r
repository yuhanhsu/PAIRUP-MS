##########################################################################################
### Generate ROC curve and AUC for each signal x met set annotation matrix
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ROCR)

# read in parameters from input argument file name
source(commandArgs(trailingOnly=T)[1],echo=T)

# ----------------------------------------------------------------------------------------
# number of input files (e.g. corresponding to different annotation matrices)
nFiles <- length(dataFiles)

# rainbow palette for plotting
col.rainbow <- rainbow(nFiles)

# function to calculate p-value for Spearman correlation
spearmanp <- function(r,n) {
	t <- r * sqrt((n-2)/((1-r)*(1+r)))
	return(2*pt(-abs(t), n-2))
}

# iterate through each input data file to plot ROC curve and calculate AUC
pdf(outFile,height=6,width=6)
legendLabs <- rep("",nFiles)
for (i in 1:nFiles) {
	df <- read.table(dataFiles[i],header=T,sep="\t")

	pValue <- spearmanp(df$Score,nMCs[i])
	pred <- prediction(pValue,df$Label)
	roc <- performance(pred,measure="tpr",x.measure="fpr")
	auc <- performance(pred,measure="auc")
	aucValue <- unlist(slot(auc,"y.values"))
	legendLabs[i] <- paste("File=",i,", AUC=",format(aucValue,digits=3),sep="")
	
	if (i == 1) {
		plot(roc,col=col.rainbow[i])
	} else {
		plot(roc,add=T,col=col.rainbow[i])
	}
}

# plot null ROC curve using permuted scores
pValue <- sample(pValue)
pred <- prediction(pValue,df$Label)
roc <- performance(pred,measure="tpr",x.measure="fpr")
auc <- performance(pred,measure="auc")
aucValue <- unlist(slot(auc,"y.values"))

plot(roc,add=T,col="gray")
legendLabs[nFiles+1] <- paste("Null, AUC=",format(aucValue,digits=3),sep="")
legendColors <- c(col.rainbow,"gray")

# add legend with AUC values
legend(0.5,0.5,legendLabs,lty=c(1,1),lwd=c(1,1),col=legendColors,box.lwd = 0,box.col = "white")

dev.off()

