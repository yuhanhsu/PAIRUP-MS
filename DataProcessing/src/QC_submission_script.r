##########################################################################################
### QC pipeline for metabolomics data
###
### authors: Claire Churchhouse, Yu-Han Hsu
##########################################################################################

rm(list=ls())

# load libraries
library(moments)

# read in the parameter file
file.specifyparams <- commandArgs(trailingOnly=T)[1]
source(file=file.specifyparams, echo=T)

# create a function to write comments to a log file
output.log <- function(somestring) {
	sink(paste(dir.out, "/logfile.txt", sep=""), append=T)
	cat(somestring)
	cat("\n")
	sink()
}

# create output directory if it doesn't exist
if(!file.exists(dir.out)) { 
	dir.create(dir.out)
}

# document parameters used in log file
output.log("QC pipeline parameters:")
output.log(" ")
output.log(paste("file.raw:", file.raw))
output.log(paste("met.type:", met.type))

output.log(paste("met.column:", met.column))
output.log(paste("IS.name:", IS.name))
output.log(paste("PP.label", PP.label))
output.log(paste("biorec.label:",biorec.label))

output.log(paste("sp.scale:", sp.scale))
output.log(paste("do.extrap:", do.extrap))
output.log(paste("sp.PPnum:", sp.PPnum))
output.log(paste("sp.window:", sp.window))
output.log(paste("tol.spline:", tol.spline))

output.log(paste("plot.it:", plot.it))
output.log(paste("plot.allsplines:", plot.allsplines))
output.log(paste("n.splineplot:", n.splineplot))

output.log(paste("cv.PP:", cv.PP))
output.log(paste("sd.point:", sd.point))
output.log(paste("sd.window:", sd.window))
output.log(paste("file.sampleRunWindow:", file.sampleRunWindow))

output.log(paste("plot.final:", plot.final))

output.log(" ")
output.log(paste("The output directory is:", dir.out))
output.log("--------------------------------------------------")
output.log(" ")

# ----------------------------------------------------------------------------------------
# Step 1: Initial cleanup and normalization to internal standards and control samples

source(file=paste(dir.src, "/QC_step1_normalization.r", sep=""), echo=T)

# ----------------------------------------------------------------------------------------
# Step 2: Post-normalizaiton filtering

source(file=paste(dir.src, "/QC_step2_postnormFiltering.r", sep=""), echo=T)

# ----------------------------------------------------------------------------------------
# look at properties of filtered metabolites

source(file=paste(dir.src, "/QC_filteredInfo.r", sep=""), echo=T)


output.log("QC COMPLETED")
