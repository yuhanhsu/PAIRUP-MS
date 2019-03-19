# PAIRUP-MS

PAIRUP-MS (Pathway Analysis and Imputation to Relate Unknowns in Mass Spectrometry-based metabolite data) is a suite of computational methods for analyzing unknown signals in mass spectrometry (MS)-based untargeted metabolomics datasets. PAIRUP-MS consists of three main components: (1) a data processing pipeline for cleaning, processing, and preparing the metabolite data for downstream analyses, (2) an imputation-based approach for pairing up unknown signals across two datasets, which would enable meta-analysis of matched signals across studies, and (3) a pathway annotation and enrichment analysis framework that links the unknown signals to plausible biological functions, without needing to confirm their chemical identities.

# Documentation

Most current documentation with toy example: [PAIRUP-MS_v1.2_Documentation.pdf](PAIRUP-MS_v1.2_Documentation.pdf)

# Software and Hardware Requirements

All scripts were tested using R v3.2 or Python v2.7. The R scripts make use of the following R packages (and their dependencies): moments (v0.14), lawstat (v3.0), plyr (v1.8.4), mice (v2.25), abind (v1.4-3), ggplot2 (v2.2.1), and ROCR (v1.0-7). Some Python scripts require scipy (v0.16.0).

Memory and run time will depend on the size of your datasets. For reference, for our largest dataset (containing ~580 samples and ~15,000 metabolite signals), DataProcessing scripts took a total of ~1 day to run, ImputationMatching took < 12 hours, and PathwayAnalysis took ~2 days (with “ConstructAnnotationMatrix.py” taking most of the time). Maximum memory requirement was < 14 GB.

# Citation

Hsu, YH et al. PAIRUP-MS: Pathway analysis and imputation to relate unknowns in profiles from mass spectrometry-based metabolite data. PLOS Comput. Biol. (2019). [link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006734)

# Troubleshooting

Please contact Yu-Han Hsu (yuhanhsu@broadinstitute.org) with a brief description of your problem.
