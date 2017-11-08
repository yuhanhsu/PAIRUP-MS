#!/bin/bash

##########################################################################################
### Shell script to test run all ImputationMatching scripts
##########################################################################################

echo "Imputing signals across datasets..."
Rscript ../src/ImputeSignals.r example_ImputeSignals_param.r

echo ""
echo "Matching signals using imputation"
Rscript ../src/MatchSignalsByImp.r example_MatchSignalsByImp_noAdduct_param.r
Rscript ../src/MatchSignalsByImp.r example_MatchSignalsByImp_adduct_param.r

echo ""
echo "Matching signals using RT"
Rscript ../src/MatchSignalsByPredRt.r example_MatchSignalsByPredRt_noAdduct_param.r
Rscript ../src/MatchSignalsByPredRt.r example_MatchSignalsByPredRt_adduct_param.r

echo ""
echo "Generating COMBINED matches..."
# (can only be run if corresponding "noAdduct" and "adduct" matches were genereated)
python ../src/GetCombinedMatches.py example_GetCombinedMatches.cfg

echo ""
echo "Generating UNIQUE and RECIPROCAL matches..."
python ../src/GetUniqueReciprocalMatches.py example_GetUniqueReciprocalMatches.cfg

echo ""
echo "Calculating statistics for matched shared knowns..."
Rscript ../src/GetMatchedKnownStats.r example_GetMatchedKnownStats_param.r

echo "TEST RUN COMPLETED"

