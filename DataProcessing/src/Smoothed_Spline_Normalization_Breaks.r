##########################################################################################
### Function to perform spline normalization (called by "QC_step1_normalization.r")
###
### author: Claire Churchhouse
##########################################################################################

Plot_Spline_Smoothing <- function(Metab.Set, dir_out, Y, spln.scale, spl.PPnum, spl.window, tol.x, plot.allsplines, bk.jump) {
	
	output.log("Running the spline smoothing function.")
	
	h <- grep(biorec.label, colnames(Y))
	k <- grep(PP.label, colnames(Y))
	CLR <- rep(1, ncol(Y))
	CLR[h] <- 3; CLR[k] <- 2
	PCH <- rep(1, ncol(Y))
	PCH[h] <- 6; PCH[k] <- 15
	
	# instead of user defined spline endpoints, consider windows of spl.window points,
	# and fit splines on to the spl.PPnum nearest PP points to infer the value at that window.
	n.fullwin <- floor(ncol(Y)/spl.window)
	
	# spline drawing palette
	clr.spline <- rainbow(n.fullwin)[sample(n.fullwin, n.fullwin, replace=F)]
	
	firstbit <- rep(1:n.fullwin, each=spl.window)
	index.x1 <- c(firstbit, rep(n.fullwin, ncol(Y)-length(firstbit)))
	# create a matrix of distances between each PP sample and the midpoint of these windows:
	mdpts <- rep(0, n.fullwin)
	for (k in 1:n.fullwin) {
		mdpts[k] <- floor(median(which(index.x1 == k))) # identify the midpoints
	}
	# find the distance from each PP sample to each midpoint
	d.mdptPP <- matrix(0, nrow=n.fullwin, ncol=length(grep(PP.label, colnames(Y))))
	PP.nearest <- matrix(0, nrow=n.fullwin, ncol=spl.PPnum)
	for (k in 1:n.fullwin) {
		d.mdptPP[k,] <- abs(grep(PP.label, colnames(Y)) - mdpts[k])
		PP.nearest[k,] <- sort(order(d.mdptPP[k,])[1:spl.PPnum])
	}
	
	if (plot.allsplines) { # if you want to plot the splines for all metabolites
		output.log("Plotting spline fits for all metabolites.")
		plot.list <- Metab.Set
	} else { # if you want to randomly select n metabolites to plot.
		output.log("Plotting spline fits for a random sample of metabolites.")
		plot.list <- sort(sample( Metab.Set, n.splineplot)) 
	}

	# if you are going to use the fwd_vs_bckw predictor determined breakpoints, then read these in
	if (bk.jump==F) {
		output.log("Using forward versus backward predictor differences to detect breakpoints.")
		nf <- max(count.fields(nearPP.file))
		list.PP <- as.list(as.data.frame(t(read.table(nearPP.file, as.is=T, skip=1, fill=T, col.names=1:nf,
			colClasses=c(rep("NULL", 4), rep("character", nf-4)))))) 
		List.PP <- lapply(list.PP, function(x) x[!is.na(x)])
		names.PP <- as.list(as.data.frame(t(read.table(nearPP.file, as.is=T, skip=1, fill=T, col.names=1:nf,
			colClasses=c("NULL", "character", rep("NULL", nf-2)))))) 
		names(List.PP) <- as.character(unlist(names.PP))
	}
	
	output.log(paste("Normalizing", length(Metab.Set), "metabolites."))
	OUT <- NULL
	output.log("Spline fitting underway...")
	for (i in Metab.Set) { # i is the row number of Y

		y.m <- log10(Y[i,]) 
		PP <- grep(PP.label, names(y.m))
		
		# calculate the mean change in PP values between consecutive PP (account for the fact there are missing PPs)
		if (bk.jump==TRUE) {
			# use 3sd jumps in PP as breakpoint indicator
			dummy <- y.m[PP]
			dummy <- dummy[!is.na(dummy)]
			if (length(dummy)==0) { 
				
				# if all PP are NA, skip this metabolite
				output.log(paste("All PP are NA: Skipping spline normalization of metabolite #",i,colnames(Y)[i]))
				OUT <- rbind(OUT, NA)
				rownames(OUT)[nrow(OUT)] <- rownames(Y)[i]
				write.table(paste(i, rownames(Y)[i], "NA_PP", sep="\t"),
					file=paste(dir_out,"/Filtering_metabolites.txt",sep=""), row=F, col=F, quote=F, append=T)
				next
			} 
			PP.jump <- abs(dummy[2:length(dummy)] - dummy[1:(length(dummy)-1)])
			PP.jump.norm <- (PP.jump - mean(PP.jump))/sd(PP.jump)
			# if normalized jump between consecutive PP samples is more than 3,
			# then insist there is a breakpoint at this PP sample
			breakPPs <- which(PP.jump.norm > 3)
			Nbreaks <- (length(breakPPs)-1)
			# do not allow consecutive points to be breakpoints, in case the first one is an outlier.
			# remove consecutive breaks
			if (Nbreaks > 0) {
				w <- 1
				while(w <= Nbreaks) {
					if(breakPPs[w+1] == (1+breakPPs[w])) {
						breakPPs <- breakPPs[-(w+1)]
						w <- w + 1
						Nbreaks <- (length(breakPPs)-1)
					} else { w <- w + 1 }
				}
			}
			break.here <- names(breakPPs)			
		} else {
			# use the fwd_vs_bckw determined breakpoints
			wn <- which(names(List.PP) ==  gsub("/", "_", gsub(" ", "", rownames(Y)[i])))
			break.all <- as.character(unlist(List.PP[wn]))
			break.here <-  break.all[which( break.all!="")]
		}
			
		break.indx <- NULL
		for (q in break.here) {
			break.indx <- c(break.indx, which(q == colnames(Y)))
		}
		if (length(break.here) >= 1) {
			for (b in 1:length(break.here)) {
				write.table(paste(i, rownames(Y)[i], break.indx[b], break.here[b], "Break", sep="\t"),
					file=paste(dir_out, "/Filtering_points.txt", sep=""), row=F, col=F, quote=F, append=T)
			}
		}
		# specify these breakpoints by adjusting PP.nearest for this metabolite i
		PP.nearest.i <- PP.nearest
		# for each breakpoint, add NAs to PP.nearest to implement the break so they don't get used in spline fitting
		if (length(break.here) >= 1) {
			# for each breakpoint
			for (q in break.here) {
				# does each window (i.e for each spline)
				for (k in 1:nrow(PP.nearest.i)) {
					these.PP <- colnames(Y)[grep(PP.label, colnames(Y))[PP.nearest.i[k,]]]
					if (is.element(q, these.PP)) {
						this.point <- which(q == these.PP) 
						if (this.point > (ncol(PP.nearest.i)/2)) { #is the breakpoint to the left or right of the window center
							PP.nearest.i[k,(this.point:ncol(PP.nearest.i))] <- NA # NAs after the break point
						} else {
							PP.nearest.i[k, 1:this.point] <- NA # NAs before the break point
						}
					}
				} # end of loop over spline window
			} # end of loop over break point
		}	
		
		# exclude outlier PP samples from the spline fitting
		y.normedPP <- (y.m[PP] - mean(y.m[PP], na.rm=T)) / sd(y.m[PP], na.rm=T)
		y.holder <- y.m[PP]
		sum.outlierPP <- sum(abs(y.normedPP) > 3, na.rm=T)
		if (sum.outlierPP > 0) { 
			outPP <- which((abs(y.normedPP) > 3)=="TRUE")
			for (b in 1:length(outPP)) {	
				write.table(paste(i, rownames(Y)[i], PP[b], names(y.m[PP])[b], "Outlier_PP", sep="\t"),
					file=paste(dir_out, "/Filtering_points.txt", sep=""), row=F, col=F, quote=F, append=T)
			}
		}
		y.holder[which(abs(y.normedPP) > 3)] <- NA
		y.m[PP] <- y.holder
	
		# cannot pass NA to smooth.spline so remove them from PP
		PP <- PP[!is.na(y.m[PP])]
		
		# fit each spline in turn
		spl.out.y <- NULL
		for (j in 1:n.fullwin) {	
			
			range.x1 <- which(index.x1 == j) # a window of length spl.window samples
			PP1 <- (grep(PP.label, names(y.m)))[PP.nearest.i[j,]] # the spl.PPnum closest PP samples to the midpoint of the rnage.x1 window
			PP1 <- PP1[!is.na(y.m[PP1])] # cannot pass NA to smooth.spline

			nx <- sum(!duplicated(round((PP1 - mean(PP1))/tol.x)))
			df.1 <- nx*spln.scale
			df.1 <- min(nx, df.1) # df must be between 1 and number of unique x
				
			if (nx >= 4) { # need at least 4 x values to fit a smoothed spline
				print(paste("Metabolite", i, "spline", j, "df.1=", df.1))
				spl.out1 <- smooth.spline(x=PP1, y=y.m[PP1], df=df.1, tol=tol.x)
				spl.pred1 <- predict(spl.out1, x=range.x1)
				
				# if opting to not extrapolate outside of PP endpoints, set these to NA
				if (!do.extrap) {
					extrap <- setdiff(range.x1, PP1[1]:PP1[length(PP1)])
					here.extrap <- NULL
					for(m in extrap){
						here.extrap <- c(here.extrap, which(range.x1 == m))
					}
					spl.pred1$y[here.extrap] <- NA
				}
				spl.pred1.y <- spl.pred1$y 
				assign(paste("spl.pred.",j,sep=""), spl.pred1.y)
				assign(paste("range.x.",j,sep=""), range.x1)
				rm(spl.pred1)
			} else {
				spl.pred1.y <- rep(NA, length(range.x1))
				assign(paste("spl.pred.", j, sep=""), spl.pred1.y)
				assign(paste("range.x.", j, sep=""), range.x1)
			}
			spl.out.y <- c(spl.out.y, spl.pred1.y) # NB. this doesn't get used in plotting
			rm(PP1, range.x1, df.1, spl.pred1.y); # rm(rangePP1)
		}

		# check that the fitted values are the correct length:
		if (length(spl.out.y) != ncol(Y)) {
			stop(paste("Length of fitted spline not equal to total number of samples!"))
		}

		# normalize to the mean of the fitted spline at the PP samples
		mn <- mean(spl.out.y[PP] , na.rm=T)
		fct <- mn / spl.out.y
		y.splinePP <- fct * y.m	
	
		# if you opt to plot the results	
		notallNA <- (sum(is.na(y.splinePP)) != length(y.splinePP))
		if (plot.it & is.element(i, plot.list) & notallNA) {
		
			# create output directory if it doesn't exist
			dir.spline.out <- paste(dir_out, "/spline_plots", sep="")
			if (!file.exists(dir.spline.out)) { 
				dir.create(dir.spline.out)
				output.log(paste("Creating the spline plot directory", dir.spline.out))
			}
			
			title.text <- paste(met.type, ": Row", i, "Compound", rownames(Y)[i])
			pdf(file=paste(dir.spline.out, "/", met.type, "_Spline_Row", i, "_",
				gsub("/", "_", rownames(Y)[i]), ".pdf", sep=""), height=9.3, width=16)
			par(mfrow=c(2,1))
			
			# plot the IS norm values, before any PP spline normalization 
			plot(1:ncol(Y), log10(Y[i,]), main="Internal standard normalized",
				xlab="run order", ylab="log_10-abundance", col=CLR, pch=PCH)
			# if a PP spline was fitted, draw it in
			for (j in 1:n.fullwin) {
				if (exists(paste("spl.pred.",j,sep=""))) {
					lines(x=get(paste("range.x.", j, sep="")), y=get(paste("spl.pred.", j, sep="")),
						lty=1, col=clr.spline[j], lwd=3)
				}
			}
			# draw where the breakpoints have been infered	
			abline(v=break.indx, lty=2, col=4)	
				
			# plot the spline normalized values	
			if (sum(!is.na(y.splinePP)) > 0) {
				plot(1:length(y.splinePP) , y.splinePP, main="Smoothed spline normalized",
					xlab="run order", ylab="log_10-abundance", col=CLR, pch=PCH)
			}
			# draw where the breakpoints have been infered	
			abline(v=break.indx, lty=2, col=4)	
			mtext(text=title.text, outer=TRUE, side=3, cex=1.3)
			dev.off()
		}
	
		OUT <- rbind(OUT, y.splinePP)
		rownames(OUT)[nrow(OUT)] <- rownames(Y)[i]
		rm(y.splinePP, y.m, spl.out.y, mn, fct, list=ls(pattern="^spl.pred.")) 
	}
	output.log("Spline fitting complete.")

	return(OUT)	

}


