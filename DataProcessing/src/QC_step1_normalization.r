##########################################################################################
### QC Step 1: Initial cleanup and normalization to internal standards and PP samples
###
### authors: Claire Churchhouse, Yu-Han Hsu 
##########################################################################################

output.log("Opening script QC_step1_normalization.r")

# read in the un-normalized abundance data 
X.raw  <- read.table(file.raw, header=T, sep="\t", as.is=T)
output.log(paste("Read in file:", file.raw))
X.deets <- X.raw[,1:4]

# create metabolite name labels for untargeted peaks
met.col <- which(colnames(X.deets)==met.column)
output.log(paste("Reading metabolite names from column:", met.col))
here <- (X.deets[, met.col]=="")
num.untarg.metabs <- sum(here)
num.targ.metabs <- nrow(X.raw) - sum(here)
rownames(X.raw)[!here] <- X.deets[!here, met.col]
rownames(X.raw)[here] <- paste(met.type, X.deets[here,"Compound"], sep="_")
X.deets[here, met.col] <- rownames(X.raw)[here]
write.table(X.deets, file=paste(dir.out, "/Metabolite_info.txt", sep=""), row=F, col=T, sep="\t", quote=F)

# replace any zero values with NAs (missing values)
is.zero <- (X.raw[,5:ncol(X.raw)]==0)
X.raw[,5:ncol(X.raw)][is.zero] <- NA
output.log(paste(sum(is.zero, na.rm=T), "zeros were replaced with NAs."))

# write to file the summary of each metabolite
for(i in 1:nrow(X.raw)){
	if(i == 1){
		write.table(t(c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","NAs")),
			file=paste(dir.out, "/Metabolite_raw_summaries.txt", sep=""), quote=F, sep="\t", row=F, col=F)
	}
	out.summary <- as.vector(summary(unlist(log10(X.raw[i,5:ncol(X.raw)]))))
	write.table(t(out.summary), file=paste(dir.out, "/Metabolite_raw_summaries.txt", sep=""),
		append=TRUE, quote=F, sep="\t", row=TRUE, col=F)
}

# count number of experimental samples and PP samples
X.raw <- X.raw[,5:ncol(X.raw)]
sig <- grep(PP.label, colnames(X.raw))
expm <- grep("X", colnames(X.raw) )
X.raw <- as.matrix(X.raw)
output.log(paste(length(expm), "experimental samples and", length(sig), "pooled plasma controls."))

# filter peaks with 50% or more missing samples in any quarter of the run
wdth <- floor(ncol(X.raw)/4)
A <- seq(1, ncol(X.raw), by=wdth)[1:4]
B <- A + (wdth-1)
B[4] <- ncol(X.raw)
fail.perc.na <- matrix(NA, nrow=nrow(X.raw), ncol=4)
for(i in 1:4){
	perc.na <- apply(is.na(X.raw[,A[i]:B[i]]), 1, mean)
	fail.perc.na[,i] <- (perc.na > 0.5)
	rm(perc.na )
}
remove.perc.na <- (fail.perc.na[,1] | fail.perc.na[,2] | fail.perc.na[,3] | fail.perc.na[,4])
miss50 <- rownames(X.raw)[remove.perc.na]
miss50.here <- which(remove.perc.na)
output.log(paste(length(miss50), "peaks/features removed for having at least 50% missing samples in any quarter of the run."))
if(length(miss50.here) > 0 ){
	miss.out <- data.frame(miss50.here, X.deets[miss50.here,], "50miss")
	write.table(miss.out, file=paste(dir.out, "/Filtering_metabolites.txt", sep=""),
		row=F, col=F, quote=F, append=TRUE, sep="\t")
	X.raw <- X.raw[-miss50.here,]
}

# normalize all samples (experimental and PP) to the mean IS value
output.log("Normalizing samples to IS.")
met.IS <- which(rownames(X.raw) == IS.name)
mn.IS <- mean(X.raw[met.IS,], na.rm=T) 
ft.IS <- (mn.IS / X.raw[met.IS,]) # normalize every metabolite of each sample by this factor
X.IS <- matrix(NA, nrow=nrow(X.raw), ncol=ncol(X.raw))
for(m in 1:nrow(X.raw)) {
	X.IS[m,] <- X.raw[m,] * ft.IS
}
colnames(X.IS) <- colnames(X.raw)
rownames(X.IS) <- rownames(X.raw)

# filter any samples whose IS value is more than 4 sd from mean
sd.IS <- sd(X.raw[met.IS,], na.rm=T) 
outlier.IS <- (X.raw[met.IS,] > (mn.IS + 4*sd.IS)) | (X.raw[met.IS,] < (mn.IS - 4*sd.IS))
output.log(paste(sum(outlier.IS), "samples removed for being IS outliers (more than 4 sd from mean IS)."))
if(sum(outlier.IS) > 0) {
	output.log("Samples removed:")
	output.log(colnames(X.raw)[outlier.IS])
	for(mm in 1:length(outlier.IS)) {
		write.table(paste(mm, colnames(X.raw)[outlier.IS[mm]], "Outlier_IS", sep="\t"),
			file=paste(dir.out, "/Filtering_samples.txt", sep=""), row=F, col=F, quote=F, append=T)
	}
}

X.IS <- X.IS[, !outlier.IS] # keep samples that are not outliers
sig <- grep(PP.label, colnames(X.IS))
expm <- grep("X", colnames(X.IS))

# filter peaks for which the PP samples are all lower or higher than the experimental samples
X.ranked <- apply(X.IS, 1, rank, na.last="keep")
X.ranked <- t(X.ranked)
rank_outlier <- function(v, cntl, expmtl){
	x <- v[cntl] # ranks of PP controls
	y <- v[expmtl] # ranks of experimental samples
	return((max(x, na.rm=TRUE) < min(y, na.rm=TRUE)) | (min(x, na.rm=TRUE) > max(y, na.rm=TRUE)))
	# if the ranks do not overlap at all then rank_outlier is TRUE
}

PP.strat <- apply(X.ranked, 1, rank_outlier, c(sig, expm))
if(sum(PP.strat) > 0) { 
	PP.strat.here <- NULL
	for(mm in which(PP.strat)) { 
		PP.strat.here <- c(PP.strat.here, which(X.deets[,met.col] == rownames(X.ranked)[mm]))			
	} 
	PP.strat.out <- data.frame(PP.strat.here, X.deets[PP.strat.here,], "PP.strat")
	write.table(PP.strat.out, file=paste(dir.out, "/Filtering_metabolites.txt", sep=""),
		row=F, col=F, quote=F, append=TRUE, sep="\t")
	X.IS <- X.IS[-which(PP.strat),]
}
output.log(paste("Removed", sum(PP.strat), "peaks for having PPs all smaller or all greater than the experimental samples."))


# ----------------------------------------------------------------------------------------
# Calculate the number of breakpoints across the run using the difference between 
# forward and backward linear predictors as a an indicator of a breakpoint.

# y is the vector of metabolite values; w is the length of the window over which a linear predictor will be built;
# d is the difference threshold between forward and backward predictors at which we deduce a breakpoint; 
# n is the number of samples between samples at which we will apply the predictors

dir.create(path=paste(dir.out, "/BreakpointDetection", sep="") , mode="0755")

fwd_bkwd_diff <- function(y, x, w, n, lv) {
	out <- NULL
	Midpt <- NULL
	for(i in seq(2*w, (length(y)-2*w), by=n)) {
		a <- i:(i+w-1)
		xa <- x[a]
		midpt <- x[max(a)+1]
		Midpt <- c(Midpt, midpt)
		metab.a <- y[a]
		b <- (max(a)+2):(max(a)+2+w-1)
		xb <- x[b]
		metab.b <- y[b]
		if(sum(!is.na(y[a])) > (w/2) & sum(!is.na(y[b])) > (w/2)) {
			fwd.lm <- lm(metab.a~xa, na.action="na.exclude")
			fwd.pred <- predict(object=fwd.lm, newdata=data.frame(xa=midpt), interval="prediction", level=lv)
			bckw.lm <- lm(metab.b~xb, na.action="na.exclude")
			bckw.pred <- predict(object=bckw.lm, newdata=data.frame(xb=midpt), interval="prediction", level=lv)	
			pred.lims <- rbind(c(fwd.pred[,"lwr"], fwd.pred[,"upr"]), c(bckw.pred[,"lwr"], bckw.pred[,"upr"]))
			smallest.lwr <- which(rank(c(fwd.pred[,"lwr"], bckw.pred[,"lwr"]))==1)
			overlap <- pred.lims[smallest.lwr, 2] >= pred.lims[3-smallest.lwr, 1]
		} else{ overlap <- NA }
		out <- c(out, overlap)
	} 
	return(data.frame(Midpt, out))
}

PP.sam <- grep(PP.label, colnames(X.IS))
bio.sam <- grep(biorec.label, colnames(X.IS))
exp.sam <-  sort(setdiff(1:ncol(X.IS), c(PP.sam, bio.sam)), decreasing=F)
x <- exp.sam

CLR <- rep(1, ncol(X.IS))
CLR[PP.sam] <- 2
CLR[bio.sam] <- 3

bk.file <- paste(dir.out, "/BreakpointDetection/Breakpoint_locations.txt", sep="")
nearPP.file <- paste(dir.out, "/BreakpointDetection/NearPP_locations.txt", sep="")

w <- 10; n <- 1; lv <- 0.95
rem.bkpt <- NULL
output.log("Beginning breakpoint detection...")
for(i in (1:nrow(X.IS))[-met.IS] ) {
	name.edited <- gsub("/", "_", gsub(" ", "", rownames(X.IS)[i]))
	print(paste("Breakpoint detection", name.edited))
	y <- log10(X.IS[i, exp.sam])
	bp <- fwd_bkwd_diff(y[!is.na(y)], x[!is.na(y)], w, n, lv) # calculate difference in fwd and bkwd predictors at each point
	these <- which(!bp$out)
	
	if(length(these) > 0) {
		x.these <- bp$Midpt[these]
	
		# determine how many of these are within 20 samples of each other
		j <- 2; N <- 1; indx <- NULL; indx[1] <- 1
		while(j <= length(x.these)) { # group samples that exceed ylm by their run order position if they are within 20 samples of each other
			if(abs(x.these[j]-x.these[j-1]) < 20) { indx[j] <- N } else { N <- N + 1; indx[j] <- N }
			j <- j + 1
		}
		num_breakpts <- length(unique(indx)) # calculate number of breakpoints

		# for each within-20 group of points, take the average of the group to declare as the breakpoint
		j <- 1; bkp <- rep(NA, length(unique(indx)))
		for(j in unique(indx)) {
			bkp[j] <- floor(mean(x.these[indx==j]))
		}

		# write the location of breakpoints to file 
		if(!file.exists(bk.file)) {
			write.table(t(c("row", "metabolite", "w", "level", "breakpoints")),
				file=bk.file , sep="\t", row=F, col=F, quote=F)
		}
		write.table(t(c(i, name.edited, w, lv, bkp)),
			file=bk.file , sep="\t", row=F, col=F, quote=F, append=T)
	
		# write the name of the nearest PP sample to each breakpoint to file
		near.PP <- NULL
		for(j in 1:length(bkp)) {
			near.PP <- c(near.PP, colnames(X.IS)[PP.sam[order(abs(PP.sam - bkp[j]), decreasing=F)[1]]])
		}
		if(!file.exists(nearPP.file)) {
			write.table(t(c("row", "metabolite", "w", "level", "breakpoints")),
				file=nearPP.file , sep="\t", row=F, col=F, quote=F)
		}
		write.table(t(c(i, name.edited, w, lv, near.PP)), file=nearPP.file, sep="\t", row=F, col=F, quote=F, append=T)
			
		# if 3 or more breakpoints are detected, then record the index to be removed from the dataset
		if(length(bkp) >= 3) {
			rem.bkpt <- c(rem.bkpt, i)
		}
		rm(num_breakpts, x.these, bkp, N)
	} else {
		 # if no breakpoints are detected
		 # write the location of breakpoints as "NA" to file 
		if(!file.exists(bk.file)) {
			write.table(t(c("row", "metabolite", "w", "level", "breakpoints")),
				file=bk.file , sep="\t", row=F, col=F, quote=F)
		}
		write.table(t(c(i, name.edited, w, lv, "NA")), file=bk.file , sep="\t", row=F, col=F, quote=F, append=T)
	
		# write "NA" to file
		if(!file.exists(nearPP.file)) {
		write.table(t(c("row", "metabolite", "w", "level", "breakpoints")),
			file=nearPP.file , sep="\t", row=F, col=F, quote=F)
		}
		write.table(t(c(i, name.edited, w, lv, "NA")), file=nearPP.file , sep="\t", row=F, col=F, quote=F, append=T)	
	}
	rm(y, bp, these)
}
output.log("Completed breakpoint detection.")

# remove the peaks with 3 or more bkpts
output.log(paste("Removed", length(rem.bkpt), "peaks for having 3 or more breakpoints."))
if(length(rem.bkpt) > 0) {
	rem.bkpt.here <- NULL
	for(k in rem.bkpt) {
		rem.bkpt.here <- c(rem.bkpt.here, which( X.deets[,met.col] == rownames(X.IS)[k]))
	}
	bkpt.out <- data.frame(rem.bkpt.here, X.deets[rem.bkpt.here,], "Breaks")
	write.table(bkpt.out, file=paste(dir.out, "/Filtering_metabolites.txt", sep=""),
		row=F, col=F, quote=F, append=T, sep="\t")
	X.IS <- X.IS[-rem.bkpt,] 
}

# look at frequency of breakpoints across the run order
nf <- max(count.fields(bk.file))
list.bkpt <- as.list(as.data.frame(t(read.table(bk.file, as.is=T, skip=1, fill=T, col.names=1:nf, 
	colClasses=c(rep("NULL", 4), rep("integer", nf-4)))))) 
List.bkpt <- lapply(list.bkpt, function(x) x[!is.na(x)])
names.bkpt <- as.list(as.data.frame(t(read.table(bk.file, as.is=T, skip=1, fill=T, col.names=1:nf,
	colClasses=c("NULL", "character", rep("NULL", nf-2)))))) 
names(List.bkpt) <- as.character(unlist(names.bkpt))
out <- NULL
for (k in 1:length(List.bkpt)) {
	out <- c(out, List.bkpt[[k]]) 
}
out <- sort(out)

if (sum(!is.na(out)) == 0) {
	# no breakpoints found
	output.log("No breakpoints were found for any peaks.")
} else {
	# look at location of breakpoints
	output.log("Plotting the histogram of breakpoint locations.")
	pdf(file=paste(dir.out, "/BreakpointDetection/Hist_breakpoint_locations.pdf", sep=""), width=6, height=6)
	hist(out, breaks=100, main="Location of breakpoints", xlab="run order")
	dev.off()

	# look at number of breakpoints within peak
	output.log("Plotting the histogram of breakpoints per peak.")
	bp.peak <- NULL
	for(k in 1:length(List.bkpt)) {
		bp.peak[k] <- length(List.bkpt[[k]])
	}
	pdf(file=paste(dir.out, "/BreakpointDetection/Hist_breakpoint_perpeak.pdf", sep=""), width=6, height=6)
	hist(bp.peak, main="Number of breakpoints per peak", xlab="Count")
	dev.off()
}

# ----------------------------------------------------------------------------------------
# Normalize to smoothed splines fitted to the PP samples, on a per peak/metabolite basis

output.log("Normalizing samples to smoothed splines fitted to PP samples.")

# read in spline normalization function
source(file=paste(dir.src,"/Smoothed_Spline_Normalization_Breaks.r",sep=""), echo=T)

# run the spline normalization
met.set <- 1:nrow(X.IS) # normalize all metabolites 
X.IS.normPP <- Plot_Spline_Smoothing(met.set, dir.out, X.IS, sp.scale, sp.PPnum, sp.window, tol.spline, 
	plot.allsplines, bk.jump=F)	
	
# save the object of IS and PP normalized abundances
ISPPfilename <- paste(dir.out, "/XISPPnorm_", met.type, ".txt", sep="")
write.table(signif(X.IS.normPP, digits=4), file=ISPPfilename, sep="\t", row=T, col=T, quote=F)
output.log(paste("Writing the data matrix of IS and PP normalized abundances to file:", ISPPfilename))
assign(paste("X.n.", met.type, sep=""), X.IS.normPP)

# output the variation of the biorec controls
output.log("Writing the variance of the biorec controls to file.")
bios <- grep(biorec.label, colnames(X.IS.normPP))
biovar <- signif(apply(X.IS.normPP[,bios], 1, var, na.rm=T), 5)
write.table(biovar, file=paste(dir.out, "/Biorec_Var.txt", sep="") , row=T, col=F, quote=F, sep="\t") 

output.log("Closing script QC_step1_normalization.r")

