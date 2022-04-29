#source all the functions so they can be accessed and used later
source("/mnt/ceph/kark6289/TCGA_analysis/trios/trio.code/trios.R")
source("/mnt/ceph/kark6289/TCGA_analysis/trios/trio.code/findTrioAll.R")
source("/mnt/ceph/kark6289/TCGA_analysis/trios/trio.code/findDups.R")
source("/mnt/ceph/kark6289/TCGA_analysis/trios/trio.code/removeNA.R")

#set the working directory
setwd("/mnt/ceph/kark6289")

#load the packages
library("splitstackshape")
library("org.Hs.eg.db")

#set the directory
setwd("/mnt/ceph/megheib/Cancer_Genomics/TCGA_data_results")
getwd()

#load the library
library(data.table)

#read in the datasets

#Gene Expression dataset
gene <- fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt")
dim(gene)

#CNA dataset
cna <- fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_CNA.txt")
dim(cna)

#Load the Methylation dataset
load("TCGA.meth.RData")

#copy to a new variable
meth <- TCGA.meth
dim(meth)


################################### matching trios ##############

#go to the findTrioAll function and run the code and get the trios dataset
final <- findTrioAll(meth, cna, gene)
print(final[1:5,])

#save the trios dataset to the path in txt format
write.table(final, file = "/mnt/ceph/kark6289/TCGA_analysis/trios/final.copy.txt", sep = "\t", row.names = FALSE,
                         col.names = TRUE, append = TRUE,quote=FALSE)

final.tmp <- cna.add.dups(final, org.Hs.eg.db)
final.res <- gene..add.dups(final.tmp, org.Hs.eg.db)

#save the data to txt file
write.table(final.res, file = "/mnt/ceph/kark6289/TCGA_analysis/trios/Trios.final.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE, append = TRUE,quote=FALSE)
