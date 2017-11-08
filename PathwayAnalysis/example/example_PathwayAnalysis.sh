#!/bin/bash

##########################################################################################
### Shell script to test run all PathwayAnalysis scripts
##########################################################################################

### PART A: GENERATING PATHWAY ANNOTATIONS

echo "Creating metabolite sets..."
python ../src/CreateMetaboliteSets.py example_CreateMetaboliteSets.cfg

echo ""
echo "Deriving metabolic components..."
Rscript ../src/DeriveMetabolicComponents.r example_DeriveMetabolicComponents_param.r

echo ""
echo "Calculating MC-metabolite set enrichment..."
python ../src/CalculateMcMetSetEnrichment.py example_CalculateMcMetSetEnrichment.cfg

echo ""
echo "Getting list of MCs to use for reconstitution..."
tail -n +2 example_data_1_MC_CPDB_2Mets_SetEnrichment_fdr20.txt | \
cut -f 1 | head -10 > example_data_1_MC_CPDB_2Mets_SetEnrichment_top10MCs.txt

echo ""
echo "Constructing signal x metabolite set annotation matrix..."
python ../src/ConstructAnnotationMatrix.py example_ConstructAnnotationMatrix.cfg

echo ""
echo "Calculating post-reconstitution met set label confidence scores..."
python ../src/CalculatePostReconEnrichment.py example_CalculatePostReconEnrichment.cfg

echo ""
echo "Generating ROC curve and AUC using annotation matrix..."
python ../src/GenerateRocCurveData.py example_GenerateRocCurveData.cfg
Rscript ../src/PlotRocCurvesWithAuc.r example_PlotRocCurvesWithAuc_param.r

# ----------------------------------------------------------------------------
### PART B: PERFORMING PATHWAY ENRICHMENT ANALYSIS

echo ""
echo "Identifying trait-associated signals + generating null signal lists..."
python ../src/IdTraitSignals_Pearson.py example_IdTraitSignals_Pearson.cfg
python ../src/IdNullTraitSignals_Pearson.py example_IdNullTraitSignals_Pearson.cfg

python ../src/IdTraitSignals_Ranksum.py example_IdTraitSignals_Ranksum.cfg
python ../src/IdNullTraitSignals_Ranksum.py example_IdNullTraitSignals_Ranksum.cfg

echo ""
echo "Converting Dataset 2 signal lists to Dataset 1 signal lists..."
python ../src/GetMatchedSignalLists.py example_GetMatchedSignalLists.cfg
python ../src/GetMatchedSignalLists.py example_GetMatchedSignalLists_sharedKnowns.cfg

echo ""
echo "Calculating nominal metabolite set enrichment p-values..."
python ../src/CalculateRanksumP.py example_CalculateRanksumP.cfg
python ../src/CalculateRanksumP.py example_CalculateRanksumP_null.cfg

python ../src/CalculateFisherP.py example_CalculateFisherP.cfg

echo ""
echo "Calculating permutation p-values..."
python ../src/CalculatePermP.py example_CalculatePermP.cfg
python ../src/CalculatePermP.py example_CalculatePermP_null.cfg

echo ""
echo "Calculating FDR + Outputting final metabolite set enrichment results..."
python ../src/CalculateFDR.py example_CalculateFDR.cfg

echo "TEST RUN COMPLETED"

