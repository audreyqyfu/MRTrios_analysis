# Title: combining and constructing TCG meth data into one data set (TCGA.meth.RData) 

# Description: Downloading a package "TCGAbiolinks", then gathering the individual data, then adding individual information to it.

# Created by Mohamed Megheib

# Date: 05-07-2021

# Last updated: 9-16-2021


#==================================================================================================================
#Library
library(TCGAbiolinks, lib="/mnt/ceph/megheib/Rpackages")


#================================================================================================================
# Read TCG meth data
#load("TCGA_BRCA.RData")
BRCAMatrix <- assay(BRCARnaseqSE)
BRCAMatrix[1:10, 1:5]

# Read one of individual information file to use it in adding Gene names to TCG
Gene_probe_code= read.delim("jhu-usc.edu_BRCA.HumanMethylation450.27.lvl-3.TCGA-A7-A4SF-01A-11D-A268-05.txt")
Gene_probe_code.1= Gene_probe_code[, -2]
rownames(Gene_probe_code.1)= Gene_probe_code.1$Composite.Element.REF
head(Gene_probe_code.1)
Gene_probe_code.2= Gene_probe_code.1[, -1]
head(Gene_probe_code.2)


#Adding gene name to TCG data 
TCGA.meth=merge(Gene_probe_code.2, BRCAMatrix, by="row.names")
dim(TCGA.meth)
TCGA.meth[1:5, 1:6]

save(TCGA.meth,file="/mnt/ceph/megheib/Cancer_Genomics/TCGA.meth.RData")

load("TCGA.meth.RData")
