##########################################################################################
### Config file for running ImputeSignals.r
##########################################################################################

### Input Files

# input file containing metabolite dataset 1 (dataset 1 samples x signals, abundance z-scores)
metDataFile1 <- "example_met_data_1.txt"

# input file containing metabolite dataset 2 (dataset 2 samples x signals, abundance z-scores)
metDataFile2 <- "example_met_data_2.txt"

# input file containing mapping of known metabolites shared by the 2 datasets
metMapUniqFile <- "example_shared_knowns.txt"


### Output Files

# output file containing imputed dataset 1 (dataset 1 samples x dataset 2 signals)
out_impFile1 <- "example_imp_data_1.txt"

# output file containing imputed dataset 2 (dataset 2 samples x dataset 1 signals)
out_impFile2 <- "example_imp_data_2.txt"

