#source all the functions so they can be accessed and used later
source("/mnt/ceph/kark6289/test_trio/trios/trio.code/trios.R")
source("/mnt/ceph/kark6289/test_trio/trios/trio.code/findTrioAll.R")
source("/mnt/ceph/kark6289/test_trio/trios/trio.code/findDups.R")
source("/mnt/ceph/kark6289/test_trio/trios/trio.code/removeNA.R")
source("/mnt/ceph/kark6289/test_trio/trios/trio.code/addDups.R")

#set the working directory
setwd("/mnt/ceph/kark6289")

#load the packages
library("splitstackshape")
library(org.Hs.eg.db)

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
final <- findTrioAll(meth, cna, gene, nStartMeth, nStartGene)
print(final[1:5,])

colnames(final) <- c("Gene name", "meth.row", "cna.row", "gene.row")
final.colnames <- colnames(final)

setwd("/mnt/ceph/kark6289/test_trio/trios")
getwd()

write.table(t(final.colnames), file = "/mnt/ceph/kark6289/test_trio/trios/pre.final.txt", sep = "\t", row.names = FALSE,
                         col.names = FALSE, append = FALSE,quote=FALSE)


#save the trios dataset to the path in txt format
write.table(final, file = "/mnt/ceph/kark6289/test_trio/trios/pre.final.txt", sep = "\t", row.names = FALSE,
                         col.names = FALSE, append = TRUE,quote=FALSE)

#find the duplicates for genes with multiple entrez ids in CNA data
final.tmp <- cna.add.dups(final, cna)

final <- final.tmp[[1]]
print(dim(final))

dup.final <- final.tmp[[2]]
print(dim(dup.final))

#find the duplicates for genes with multiple entrez ids in Gene Exp data
final.res <- gene.add.dups(final, dup.final, gene)
print(final.res[1:10,])
print(dim(final.res))

write.table(t(colnames(final.res)), file = "/mnt/ceph/kark6289/test_trio/trios/Trios.final.txt", sep = "\t", row.names = FALSE,
                         col.names = FALSE, append = FALSE,quote=FALSE)


#save the data to txt file
write.table(final.res, file = "/mnt/ceph/kark6289/test_trio/trios/Trios.final.txt", sep = "\t", row.names = FALSE,
            col.names = FALSE, append = TRUE,quote=FALSE)

