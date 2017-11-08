##########################################################################################
### Config file for running DeriveMetabolicComponents.r
##########################################################################################

### Parameters

# number of samples in metabolite data (only need to save signal scores in 1st nSamples MCs)
nSamples <- 50


### Input Files

# input file containing metabolite abundance data (sample x signal matrix)
#	tab-delimited file with headers: SampleID followed by signal names
metDataFile <- "example_met_data_1.txt"


### Output Files

# output file containing signal-signal Spearman correlatrion matrix
#	gzipped tab-delimited file with metabolite names in headers and 1st column
#	will not be generated if value set to empty string ("")
corMatrixFile <- "example_data_1_MetaboliteCorMatrix.txt.gz"

# output file containing MC statistics
#	tab-delimited file with headers: MC, StandardDeviation, ProportionVariance, CumulativeProportion
mcStatsFile <- "example_data_1_McStatistics.txt"

# output file containing metabolite signal x MC score matrix
#	tab-delimited file with MCs as headers (no header for 1st column: signal name)
mcScoreMatrixFile <- "example_data_1_McScoreMatrix.txt"

