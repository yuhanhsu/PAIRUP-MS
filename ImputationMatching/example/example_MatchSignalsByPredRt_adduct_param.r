##########################################################################################
### Config file for running MatchSignalsByPredRt.r
##########################################################################################

### Input Files and Parameters

# profiling details of dataset 1 (add to lists below if more than 1 profiling method)
metListFile1 <- c("example_met_names_1.txt") # list of all metabolites to be matched (those in example_met_data_1.txt)
metInfoFiles1 <- c("example_met_info_1.txt") # file containing m/z and RT info
methods1 <- c("C8-pos") # profiling method name
ionModes1 <- c("POS") # ionization mode: "POS" (positive) or "NEG" (negative)

# profiling details of dataset 2 (add to lists below if more than 1 profiling method)
metListFile2 <- c("example_met_names_2.txt") # list of all metabolites to be matched (those in example_met_data_2.txt)
metInfoFiles2 <- c("example_met_info_2.txt") # file containing m/z and RT info
methods2 <- c("LIPID") # profiling method name
ionModes2 <- c("POS") # ionization mode: "POS" (positive) or "NEG" (negative)

# input file containing mapping of known metabolites shared by the 2 datasets
metMapFile <- "example_shared_knowns.txt"

# parameters for matching m/z values
mzWindow <- 0.005
mzMode <- "adduct" # "noAdduct" or "adduct"
posAdducts <- c(0,1.007276,18.033823,22.989218) # list of POS adducts to consider (each value corresponds to mass of X in M+X: M, M+H, M+NH4, M+Na)
negAdducts <- c(-1.007276) # list of NEG adducts to consider (e.g. M-H adduct ion)


### Output Files

# output files containing signal pairs matched using m/z + RT:
out_match12_RtFile <- paste("example_1-2_MatchedSignals_ByPredRt_",mzMode,"MzMode.txt",sep="")
out_match21_RtFile <- paste("example_2-1_MatchedSignals_ByPredRt_",mzMode,"MzMode.txt",sep="")

