##########################################################################################
### QC Step 2: Post-normalization filtering
###
### authors: Claire Churchhouse, Yu-Han Hsu
##########################################################################################

output.log("Opening script QC_step2_postnormFiltering.r")

# remove metabolite if CV of PP is greater than cv.PP:
PPcols <- grep(PP.label, colnames(X.IS.normPP))
MDN <- apply(X.IS.normPP[,PPcols], 1, median, na.rm=T)
SD <- apply(X.IS.normPP[,PPcols], 1, sd, na.rm=T)
CV <- SD/MDN
here <- which(CV > cv.PP)
if (length(here) > 0) {
	cvPP.here <- NULL
	for (mm in here) {
		cvPP.here <- c(cvPP.here, which(X.deets[,met.col] == rownames(X.IS.normPP)[mm]))	
	}
	cvPP.out <- data.frame(cvPP.here, X.deets[cvPP.here,], "cvPP")
	write.table(cvPP.out, file=paste(dir.out, "/Filtering_metabolites.txt", sep=""),
		row=F, col=F, quote=F, append=T, sep="\t")

	output.log(paste(length(here), "metabolites removed for CV of PP samples being greater than", cv.PP))
	X.QC <- X.IS.normPP
	X.QC <- X.QC[-here,]
} else {
	X.QC <- X.IS.normPP
	output.log("No metabolites removed for CV of PP samples being greater than cv.PP") 
}


# remove points that are phenotypic outliers of sd.point SD or more
mn <- apply(X.QC, 1, mean, na.rm=T)
MN <- matrix(rep(mn, ncol(X.QC)), nrow=nrow(X.QC), ncol=ncol(X.QC), byrow=F)
stdev <- apply(X.QC, 1, sd, na.rm=T)
STDEV <- matrix(rep(stdev, ncol(X.QC)), nrow=nrow(X.QC), ncol=ncol(X.QC), byrow=F)
X.std <- abs(X.QC - MN)/STDEV
nthese <- 0
for (i in 1:nrow(X.QC)) {
	these <- which(abs(X.std[i,]) > sd.point)
	if (length(these) > 0) {
		write.table(cbind(i, rownames(X.QC)[i], these, colnames(X.QC)[these], "pheno_outlier"),
			file=paste(dir.out, "/Filtering_points.txt", sep=""), sep="\t", append=T, row=F, col=F, quote=F)
		nthese <- nthese + length(these)
		X.QC[i,these] <- NA # remove these outliers
	}
}
output.log(paste(nthese, "points removed for being phenotypic outliers more than sd.point SD."))


# remove windows with variance or mean > sd.window (across all windows of a metabolite)
tempSampLabels <- gsub("^X","", colnames(X.IS.normPP))
windowGrouping <- read.table(file.sampleRunWindow, sep="\t", header=T)
windowGroup <- windowGrouping$WindowGroup[match(tempSampLabels,windowGrouping$Sample)]

nthese <- 0
nwindows <- 0
for (i in 1:nrow(X.QC)) {
	temp <- cbind(as.numeric(X.QC[i,]),windowGroup)
	colnames(temp) = c("m","grp")

	if (sum(!is.na(temp[!is.na(temp[,2]),1])) > 0) { # only do this if there are any metabolite data
		windowVars <- aggregate(m ~ grp, temp, function(x) var(x))$m
		windowVarZs <- (windowVars - mean(windowVars,na.rm=T)) / sd(windowVars,na.rm=T)

		windowMeans <- aggregate(m ~ grp, temp, function(x) mean(x))$m
		windowMeanZs <- (windowMeans - mean(windowMeans,na.rm=T)) / sd(windowMeans,na.rm=T)

		badWindows <- which(abs(windowMeanZs) > sd.window | abs(windowVarZs) > sd.window)

		if (length(badWindows > 0)) {
			nthese = nthese + 1
			nwindows = nwindows + length(badWindows)
	
			# remove points in outlier windows
			for (w in badWindows) {
				write.table(cbind(i, rownames(X.QC)[i], w, "window_outlier"), 
					file=paste(dir.out, "/Filtering_windows.txt", sep=""),
					sep="\t", append=T, row=F, col=F, quote=F)
				X.QC[i,which(temp[,"grp"]==w)] <- NA
			}
		}
	}
}

output.log(paste(nwindows, "windows in", nthese, "metabolites removed for being window mean/variance outliers (more than sd.windowSD)"))


# remove metabolite if variance in each window is found to be different by Levene's test 
library("lawstat")
grp <- windowGroup

LEV.pval <- NULL
LEV.stat <- NULL
for (i in 1:nrow(X.QC)) {
	subQc <- X.QC[i,!is.na(X.QC[i,])]
	subGrp <- grp[!is.na(X.QC[i,])]
	
	if (length(unique(subGrp[!is.na(subGrp)])) > 1) {
		lev <- levene.test(subQc, group=subGrp)
		LEV.pval[i] <- lev$p.value
		LEV.stat[i] <- lev$statistic 
	} else {
		LEV.pval[i] <- NA
		LEV.stat[i] <- NA
	}
}
bad.var <- which(LEV.pval <= (0.05/nrow(X.QC)))

badvar.here <- NULL
for(mm in bad.var){
	badvar.here <- c(badvar.here, which( X.deets[,met.col] == rownames(X.QC)[mm]))	
}
badvar.out <- data.frame(badvar.here, X.deets[badvar.here,], "Var_diff")
write.table(badvar.out, file=paste(dir.out, "/Filtering_metabolites.txt", sep=""),
	row=F, col=F, quote=F, append=T, sep="\t")

output.log(paste(length(badvar.here), "metabolites removed for window variance differing as assessed by Levene's test"))
X.var <- X.QC
X.var <- X.var[-bad.var,]


# remove metabolites with all NAs
bad.na <- which(rowSums(!is.na(X.var)) == 0)
if (length(bad.na) > 0) {
	X.var <- X.var[-bad.na,]
}


# ----------------------------------------------------------------------------------------
# output QCed Data

# write final QCed metabolite values to file, include header and row names
# remove any blank spaces from metabolite names before saving
#colnames(X.var) <- gsub("^X", "", colnames(X.var))
#colnames(X.var) <- gsub("_.*", "", colnames(X.var))
newnames <- gsub(" ", "_", rownames(X.var))
rownames(X.var) <- newnames
file.final_output_withcontrols <- gsub(".txt", "_withcontrols.txt", file.final_output)
write.table(signif(X.var,digits=5), file=file.final_output_withcontrols, sep="\t", col=T, row=T, quote=F)
output.log(paste("Wrote data for all QCed samples to file", file.final_output_withcontrols))

# write final QCed metabolite values to file but only for the experimental samples
these <- setdiff(1:ncol(X.var), c(grep(PP.label, colnames(X.var)), grep(biorec.label, colnames(X.var))))
write.table(signif(X.var[,these], digits=5), file=file.final_output, sep="\t", col=T, row=T, quote=F)
output.log(paste("Wrote data for experimental samples to file", file.final_output))

# write final QCed metabolite names to file
write.table(rownames(X.var), file=paste(dir.out, "/Metabolites_QCpass.txt", sep=""), row=F, col=F, quote=F)


# ----------------------------------------------------------------------------------------
# calculate tuning metrics

# (1) average variance/CV of window variance/mean 
var_windowVars <- NULL
cv_windowVars <- NULL
var_windowMeans <- NULL
cv_windowMeans <- NULL
for (i in 1:nrow(X.var)) {
	temp <- cbind(as.numeric(X.var[i,]),windowGroup)
	colnames(temp) = c("m","grp")

	windowVars <- aggregate(m ~ grp, temp, function(x) var(x))$m
	var_windowVars[i] <- var(windowVars,na.rm=T)
	cv_windowVars[i] <- sd(windowVars,na.rm=T)/mean(windowVars,na.rm=T)

	windowMeans <- aggregate(m ~ grp, temp, function(x) mean(x))$m
	var_windowMeans[i] <- var(windowMeans,na.rm=T)
	cv_windowMeans[i] <- sd(windowMeans,na.rm=T)/mean(windowMeans,na.rm=T)
}

# output average window statistics
output.log("")
output.log("##### QC TUNING METRICS #####")
output.log(paste("MEDIAN VARIANCE OF WINDOW VARIANCE =",median(var_windowVars,na.rm=T)))
output.log(paste("MEDIAN CV OF WINDOW VARIANCE =",median(cv_windowVars,na.rm=T)))
output.log(paste("MEDIAN VARIANCE OF WINDOW MEAN =",median(var_windowMeans,na.rm=T)))
output.log(paste("MEDIAN CV OF WINDOW MEAN =",median(cv_windowMeans,na.rm=T)))


# (2) average variance of Biorec samples

# get Biorec samples data and calculate variance for each metabolite
X.biorec <- X.var[,grep(biorec.label, colnames(X.var))]
biovar <- apply(X.biorec,1,var,na.rm=T)

output.log(paste("MEDIAN VARIANCE OF BIOREC CONTROLS =",median(biovar,na.rm=T)))


# ----------------------------------------------------------------------------------------
# optiional: plot final QC'ed values

if (plot.final == TRUE) {

	# define sample point styles
	h <- grep(biorec.label, colnames(X.var))
	k <- grep(PP.label, colnames(X.var))
	CLR <- rep(1, ncol(X.var))
	CLR[h] <- 3; CLR[k] <- 2
	PCH <- rep(1, ncol(X.var))
	PCH[h] <- 6; PCH[k] <- 15

	# create plot directory
	dir.plot.out <- paste(dir.out, "/QCed_plots", sep="")
	if (!file.exists(dir.plot.out)) {
		dir.create(dir.plot.out)
		output.log(paste("Creating the QC'ed plot directory", dir.plot.out))
	}

	for (i in 1:nrow(X.var)) {
		pdf(file=paste(dir.plot.out, "/", met.type, "_QCed_", gsub( "/", "_", rownames(X.var)[i]), ".pdf", sep=""),
			height=4.65, width=16)
		plot(1:ncol(X.var), X.var[i,], main="Final QCed", xlab="run order", ylab="log_10-abundance", col=CLR, pch=PCH)
		dev.off()
	}

	output.log("Generated final QC'ed metabolite plots")
}

output.log("Closing script QC_step2_postnormFiltering.r")
