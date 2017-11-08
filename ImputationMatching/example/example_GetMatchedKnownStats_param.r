##########################################################################################
### Config file for running GetMatchedKnownStats.r
##########################################################################################

### Input Files and Parameters

# input file containing mapping of known metabolites shared by the 2 datasets
metMapFile <- "example_shared_knowns.txt"

# input file containing metabolite dataset 1 (dataset 1 samples x signals, abundance z-scores)
metDataFile1 <- "example_met_data_1.txt"

# input file containing metabolite dataset 2 (dataset 2 samples x signals, abundance z-scores)
metDataFile2 <- "example_met_data_2.txt"

# input file containing imputed dataset 1 (dataset 1 samples x dataset 2 signals)
#       generated by ImputeSignals.r
out_impFile1 <- "example_imp_data_1.txt"

# input file containing imputed dataset 2 (dataset 2 samples x dataset 1 signals)
#       generated by ImputeSignals.r
out_impFile2 <- "example_imp_data_2.txt"

# adduct ion parameter used to generate matches
mzMode <- "noAdduct" # "noAdduct", "adduct", or "combined"

# output files containing signal pairs matched using m/z + RT
#	generated by MatchSignalsByPredRt.r
out_match12_RtFile <- paste("example_1-2_MatchedSignals_ByPredRt_",mzMode,"MzMode.txt",sep="")
out_match21_RtFile <- paste("example_2-1_MatchedSignals_ByPredRt_",mzMode,"MzMode.txt",sep="")

# input files containing signal pairs matched using m/z + imputation
#	generated by MatchSignalsByImp.r
out_match12_imp1_File <- paste("example_1-2_MatchedSignals_ByImp1_",mzMode,"MzMode.txt",sep="")
out_match21_imp1_File <- paste("example_2-1_MatchedSignals_ByImp1_",mzMode,"MzMode.txt",sep="")
out_match12_imp2_File <- paste("example_1-2_MatchedSignals_ByImp2_",mzMode,"MzMode.txt",sep="")
out_match21_imp2_File <- paste("example_2-1_MatchedSignals_ByImp2_",mzMode,"MzMode.txt",sep="")
out_match12_impAll_File <- paste("example_1-2_MatchedSignals_ByImpAll_",mzMode,"MzMode.txt",sep="")
out_match21_impAll_File <- paste("example_2-1_MatchedSignals_ByImpAll_",mzMode,"MzMode.txt",sep="")


### Output Files

# output files summarizing shared knowns matched using m/z + RT
out_knownMatch12_rtFile <- paste("example_1-2_MatchedKnowns_ByPredRt_",mzMode,"MzMode.txt",sep="")
out_knownMatch21_rtFile <- paste("example_2-1_MatchedKnowns_ByPredRt_",mzMode,"MzMode.txt",sep="")

# output files summarizing shared knowns matched using m/z + imputation
out_knownMatch12_imp1_File <- paste("example_1-2_MatchedKnowns_ByImp1_",mzMode,"MzMode.txt",sep="")
out_knownMatch21_imp1_File <- paste("example_2-1_MatchedKnowns_ByImp1_",mzMode,"MzMode.txt",sep="")
out_knownMatch12_imp2_File <- paste("example_1-2_MatchedKnowns_ByImp2_",mzMode,"MzMode.txt",sep="")
out_knownMatch21_imp2_File <- paste("example_2-1_MatchedKnowns_ByImp2_",mzMode,"MzMode.txt",sep="")
out_knownMatch12_impAll_File <- paste("example_1-2_MatchedKnowns_ByImpAll_",mzMode,"MzMode.txt",sep="")
out_knownMatch21_impAll_File <- paste("example_2-1_MatchedKnowns_ByImpAll_",mzMode,"MzMode.txt",sep="")

# plots summarizing shared knowns matched using m/z + imputation
out_knownMatch12_imp1_plotFile <- paste("example_1-2_MatchedKnowns_ByImp1_",mzMode,"MzMode_Plots.pdf",sep="")
out_knownMatch21_imp1_plotFile <- paste("example_2-1_MatchedKnowns_ByImp1_",mzMode,"MzMode_Plots.pdf",sep="")
out_knownMatch12_imp2_plotFile <- paste("example_1-2_MatchedKnowns_ByImp2_",mzMode,"MzMode_Plots.pdf",sep="")
out_knownMatch21_imp2_plotFile <- paste("example_2-1_MatchedKnowns_ByImp2_",mzMode,"MzMode_Plots.pdf",sep="")
out_knownMatch12_impAll_plotFile <- paste("example_1-2_MatchedKnowns_ByImpAll_",mzMode,"MzMode_Plots.pdf",sep="")
out_knownMatch21_impAll_plotFile <- paste("example_2-1_MatchedKnowns_ByImpAll_",mzMode,"MzMode_Plots.pdf",sep="")

