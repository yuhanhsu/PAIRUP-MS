##########################################################################################
### Config file for running MissingValueImputation.r
##########################################################################################

# input file containing metabolite data
inFile <- "example_spQC_miss0.25.txt"

# output file name
outFile <- "example_spQC_miss0.25_medianMiceImp.txt"


### MICE parameters:

# number of MICE imputation runs
nChains <- 5

# type of regression model to use ("norm" = Bayeslian linear model)
miceMethod <- "norm"

# number of iterations per run
# (we used 10-15 for our datasets)
nIters <- 10

# number of predictors (best-correlated metabolites) to use for imputing each metabolite
# (we used 50-70 for our datasets)
nPreds <- 10

