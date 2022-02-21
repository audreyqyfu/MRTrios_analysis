#set the directory
setwd("/mnt/ceph/megheib/Cancer_Genomics")
getwd()

load("ntrez_id.2.RData")
master.list <- ntrez_id.2
  
#set the directory
setwd("/mnt/ceph/megheib/Cancer_Genomics/TCGA_data_results")
getwd()

#load the library
library(data.table)
library(splitstackshape)

#read in the datasets
gene.exp <- fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt")
dim(gene.exp)


cna.TCGA <- fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_CNA.txt")
dim(cna.TCGA)


load("TCGA.meth.RData")
dim(TCGA.meth)


#counts the number of NAs in each row
#returns the rows that do not have a sum of NAs that match with the number of columns
#so basically removes the rows that have all NAs for methylation levels

gene <- gene.exp[rowSums(is.na(gene.exp[,3:1102])) != ncol(gene.exp[,3:1102]), ]
dim(gene)
gene[1:5,1:5]

cna <- cna.TCGA[rowSums(is.na(cna.TCGA[,3:1082])) != ncol(cna.TCGA[,3:1082]), ]
dim(cna)
cna[1:5,1:5]

meth <- TCGA.meth[rowSums(is.na(TCGA.meth[,5:899])) != ncol(TCGA.meth[,5:899]), ]
dim(meth)
meth[1:5,1:5]


################################### matching trios ##############
#since the data has multiple genes in one row separated by ";", we split them
split.genes <- unlist(strsplit(as.character(meth$Gene_Symbol), ';'))

#then find the unique genes from the methylation data
unique.genes <- unique(split.genes)

#to find the row in the methylation data
meth.genes <- strsplit(as.character(meth$Gene_Symbol), ';')

tmp <- NULL
findTrioAll <- function(meth.data, cna.data, gene.data, dupsList) {

    for(i in 1:unique.genes){
    
        #apply the function with the provided data
        tmp <- cbind(tmp,trios(unique.genes[i]))
  }

dim(tmp)
tmp[1:5,]

#find the rows that have NAs for each data type
meth.na.rows <- which(is.na(tmp[,2]))
cna.na.rows <- which(is.na(tmp[,3])) 
gene.na.rows <- which(is.na(tmp[,4]))

for(i in 1:length(cna.na.rows){

    entrez.cna(cna.na.rows[i])
}

for(i in 1:length(gene.na.rows){

    entrez.cna(gene.na.rows[i])
}

}
