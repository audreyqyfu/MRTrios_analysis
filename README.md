# MRTrios_analysis
This repo consists of R scripts that aim in causal network inference with Mendelian randomization for trios consisting of copy number alteration, gene expression and DNA methylation.  The R scripts here perform preprocessing, causal inference (using functions in the R package MRTrios: https://github.com/audreyqyfu/MRTrios), and downstream analyses.  The user will need to change the paths and names of the files mentioned in these scripts to suit their own needs.

The scripts should be used in the following order:

- `DataProcessing.R`: this file performs logit transformation of the methylation data and extract positive and negative ER patients from the clinical dataset.    
- `mainTrioMatch.R`: this file forms trios by integrating the CNA, methylation, and gene expression data.
- `trio.gene.type.sep.R`     
- `main.findPCs.R`: this file calculates the principal component (PC) score matrix for methylation and for gene expression data, and identifies PCs significantly associated each trio.    
- `main.analyzeTrios.R`: this file performs causal network inference for trios and their associated confounders (i.e., PCs) to infer the causal models.   
- `HumanMethProbeInfo.R`     
- `GOanalysis.R`    
