#!/bin/bash

##########################################################################################
### Shell script to test run all DataProcessing scripts
##########################################################################################

echo "Running QC..."
Rscript ../src/QC_submission_script.r example_QC_C8-pos_param.r
Rscript ../src/QC_submission_script.r example_QC_C18-neg_param.r

echo ""
echo "Filtering samples and signals based on missingness..."
Rscript ../src/MissingnessFilter.r example_MissingnessFilter_param.r

echo ""
echo "Missing value imputation..."
Rscript ../src/MissingValueImputation.r example_MissingValueImputation_param.r

echo ""
echo "Covariate adjustment + inverse normal transformation..."
Rscript ../src/AdjustCov_InvNorm.r example_AdjustCov_InvNorm_param.r

echo "TEST RUN COMPLETED"

