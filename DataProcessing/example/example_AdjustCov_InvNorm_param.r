##########################################################################################
### Config file for running AdjustCov_InvNorm.r
##########################################################################################

# input metabolite data file (after missing value imputation)
metFile <- "example_spQC_miss0.25_medianMiceImp.txt"

# input phenotype data file (optional, if needed for getting covariate data)
ptypeFile <- "example_phenotype_data.txt"

# specify any covariates to adjust for (column names in phenotype data file)
covariates <- c("Age","Sex")
#covariates <- c() # if no covariates

# output file name for adjusted + transformed metabolite data
outFile <- "example_spQC_miss0.25_medianMiceImp_adjAgeSex_invNorm.txt"

