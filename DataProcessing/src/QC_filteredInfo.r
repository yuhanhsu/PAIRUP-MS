##########################################################################################
### Look at properties of filtered metabolites
###
### authors: Claire Churchhouse, Yu-Han Hsu 
##########################################################################################

output.log("Opening script QC_filteredInfo.r")

# read in info on which metabolites have been filtered and why
filt.met <- read.table( paste( dir.out,"/Filtering_metabolites.txt", sep=""), sep="\t", header=F)
colnames(filt.met) <- c("Index", "Compound", "m.z", "RT", "Metabolite", "Filter")

# plot location of filtered metabolites in m:z versus t space
filters <- as.character(unique(filt.met[,6]))
clrs <- 1+(1:length(filters))
CLR <- rep(NA, nrow(filt.met))
for (i in 1:nrow(filt.met)) {
	CLR[i] <- which(filters == filt.met[i,6])
}
clr.pal <- c("#5FF2F5", "#F5C52A", "#8B16F7", "#F77016", "#23C436", "#2B49F0", "#E8155B")
pch.pal <- c(1:6, 8)

# annotate all filters at once
pdf(file=paste(dir.out, "/Filtered_metabolites_mz_vs_rt.pdf", sep=""), height=9.3, width=9.3)
plot(X.deets[,"RT"], X.deets[,"m.z"], col=gray(0.7), xlab="retention time (min)", ylab="mass : charge",
	main="Filtered peaks", cex=0.5)
points(filt.met[,"RT"], filt.met[,"m.z"], col=clr.pal[CLR], pch=pch.pal[CLR], cex=0.5)
legend("topright", col=clr.pal, legend=filters, pch=pch.pal)
dev.off()

# annotate one filter at a time
for (i in 1:length(filters)) {
	pdf(file=paste(dir.out, "/Filtered_metabolites_mz_vs_rt_", filters[i], ".pdf", sep=""),height=9.3, width=9.3)
	plot(X.deets[,"RT"], X.deets[,"m.z"], col=gray(0.7), xlab="retention time (min)", ylab="mass : charge",
		main="Filtered peaks", cex=0.5)
	these <- which(filt.met[,"Filter"] == filters[i])
	points(filt.met[these,"RT"], filt.met[these,"m.z"], col=clr.pal[i], pch=pch.pal[i], cex=0.5)
	legend("topright", col=clr.pal[i], legend=filters[i], pch=pch.pal[i])
	dev.off()
}

# annotate known metabolites
pdf(file=paste(dir.out, "/All_metabolites_mz_vs_rt_Known.pdf",sep=""), height=9.3, width=9.3)
plot(X.deets[,"RT"], X.deets[,"m.z"], col=gray(0.7), xlab="retention time (min)", ylab="mass : charge", main="All peaks", cex=0.5)
these <- which(!grepl(paste(met.type,"_",sep=""), X.deets[,"Metabolite"]))
points(X.deets[these,"RT"], X.deets[these,"m.z"], col=2, pch=20, cex=0.5)
legend("topright", col=c(gray(0.7),2), legend=c("unknown", "known"), pch=c(16,20))
dev.off()

# ----------------------------------------------------------------------------------------
# plot histogram of raw abundances of filtered versus not-filtered metabolites

Summ <- read.table(paste(dir.out, "/Metabolite_raw_summaries.txt", sep=""),
	fill=TRUE, header=T, as.is=T, sep="\t",row.names=NULL)
med_abun <- Summ[, "Median"]
names(med_abun) <- X.deets[,"Metabolite"]
filtered.bool <- NULL
for (i in 1:nrow(X.deets)) {
	filtered.bool[i] <- is.element(X.deets[i,"Metabolite"], filt.met[,"Metabolite"])
}

clr_ok <- rgb(0,0,1,1/4)
clr_filt <- rgb(1,0,0,1/4)
pdf(file=paste(dir.out, "/Filtered_metabolites_abundance.pdf", sep=""), height=9.3, width=9.3)
hist(med_abun[!filtered.bool], col=clr_ok , xlab="Median abundance (log10)", main="",
	xlim=c(0, ceiling(max(med_abun,na.rm=T))), breaks=30)
hist(med_abun[filtered.bool], add=TRUE, col=clr_filt, breaks=30)
legend("topright", legend=c("filtered","not filtered"), col=c(clr_filt, clr_ok), pch=15)
dev.off()


# ----------------------------------------------------------------------------------------
# tabulate the number of metabolites filtered at each step and how many are targeted vs non-targeted
count <- rep(NA, length(filters))
count_untarg <- rep(NA, length(filters))
for (i in 1:length(filters)) {
	here <- filt.met[,"Filter"]==filters[i]
	count[i] <- sum(here)
	untarg <- grep(paste(met.type,"_",sep=""), filt.met[ here, "Metabolite"]) 
	count_untarg[i] <- length(untarg)
}

write.table(paste("Total number of targeted metabolites:", num.targ.metabs),
	file=paste(dir.out, "/Filtered_metabolites_count.txt", sep=""), row=F, col=F, quote=F)

write.table(paste("Total number of untargeted metabolites:", num.untarg.metabs),
	file=paste(dir.out, "/Filtered_metabolites_count.txt", sep=""), row=F, col=F, quote=F, append=T)

write.table(paste("Total metabolites filtered:"), file=paste(dir.out, "/Filtered_metabolites_count.txt", sep=""),
	row=F, col=F, quote=F, append=T)

for (i in 1:length(filters)) {
	write.table(paste(filters[i], count[i], sep="\t"), file=paste(dir.out, "/Filtered_metabolites_count.txt", sep=""),
		row=F, col=F, quote=F, append=T)
}

write.table(paste("\nUntargeted metabolites filtered:"), file=paste(dir.out, "/Filtered_metabolites_count.txt", sep=""),
	row=F, col=F, quote=F, append=T)

for (i in 1:length(filters)) {
	write.table(paste(filters[i], count_untarg[i], sep="\t"), file=paste(dir.out, "/Filtered_metabolites_count.txt", sep=""),
		row=F, col=F, quote=F, append=T)
}


output.log("Closing script QC_filteredInfo.r")

