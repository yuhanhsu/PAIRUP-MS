##########################################################################################
### Config file for running MissingnessFilter.r
##########################################################################################

# input files containing QC'ed metabolite data (can add more files)
qcFiles <- c("QC_output_C8-pos/X_QCed_C8-pos.txt","QC_output_C18-neg/X_QCed_C18-neg.txt")

# Internal Standard metabolites
isMets <- c("Known_C8-pos_IS", # C8-pos
	    "Known_C18-pos_IS") # C18-neg

# fraction of missingness threshold to use for filteirng samples/signals
missThreshold <- 0.25

# output file name
outFile <- paste("example_spQC_miss", missThreshold, ".txt", sep="")

