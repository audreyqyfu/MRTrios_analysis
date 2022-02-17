####################################
# function to find trios
# by matching methylation probes,
# gene expression, and copy number
# alteration (CNA)
####################################
# input
#
# methyl.data: methylation data matrix; probes in rows, individuals in columns
# expn.data: gene expression data matrix; genes in rows, individuals in columns
# cna.data: cna data matrix; genes in rows, individuals in columns
# entrez.list: master list of gene names with Entrez IDs
#
# output
# trios.rows: data matrix of 4 columns: gene name, row number
#             from each input dataset
function <- findTrios (methyl.data, expn.data, cna.data, entrez.list) {
    # identify a list of unique genes
    # for methylation probes
    
    
    # identify trios by gene name (Hugo_Symbol)
    # call function trios()
    
    
    # identify additional trios by Entrez ID
    # call function entrez()
    
    
    # find additional matches using R package org.Hs.eg.db

    
    # return a data matrix for matched trios
    # in each row: gene name, row numbers from each dataset

}



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

#function to match trios (works for one gene at a time)
trios <- function (gene.name) {
  
    #assign row numbers for the list data
    #so the genes in the same list with have the same row
    g1 <- rep(seq_along(meth.genes), sapply(meth.genes, length))
    g1[1:5]
    
    #match the gene name and fine the row
    meth.row <- g1[which(unlist(meth.genes) == gene.name)]
    meth.row
    
    
    ##### error here !!!!!!!!!!!!!!
    cna.row <- which(unlist(cna$Hugo_Symbol) == gene.name)
    cna.row

    if(length(cna.row) == 0){
      cna.row = NA
    }
    
    
    
    gene.row <- which(unlist(gene$Hugo_Symbol) == gene.name)
    gene.row
    
    if(length(gene.row) == 0){
      gene.row = NA
    }
    
    # gene.row <- match(gene.name, gene$Hugo_Symbol)
    # gene.row
    
        #we save the gene name, i (which is the row number in meth data)
        #and row numbers in cna & gene data
        results <- cbind(gene.name, meth.row, cna.row, gene.row)
        
        #return the result to the function
        return(results)
      
}

#apply the function with the provided data
tmp <- trios(unique.genes[1])
dim(tmp)
tmp[1:5,]



#find the rows that have NAs for each data type
meth.na.rows <- which(is.na(tmp[,2]))
cna.na.rows <- which(is.na(tmp[,3])) 
gene.na.rows <- which(is.na(tmp[,4]))


#function to find entrez ID rows (only works for the first match)
#entrez id are duplicates
entrez <- function(cna.row,gene.row){
  
  #use the cna row but take the column of the gene row
  #column 4 has gene rows in tmp and column 3 has cna rows in tmp
  entrez.row.cna <- match(gene$Entrez_Gene_Id[as.integer(tmp[cna.row,4])], cna$Entrez_Gene_Id)
  entrez.row.cna
  
  #check if the matched row is NA or not
  #If not proceed
  if(is.na(cna$Entrez_Gene_Id[entrez.row.cna]) != TRUE){
    
  #replace the matched rows with the respective NA spots in the tmp data
  tmp[cna.row,3] <- entrez.row.cna
  
  }
  
  #use the gene row but take the column of the cna row
  entrez.row.gene <- match(cna$Entrez_Gene_Id[as.integer(tmp[gene.row,3])], gene$Entrez_Gene_Id)
  entrez.row.gene
  
  #check if the matched row is NA or not
  #If not proceed
  if(is.na(gene$Entrez_Gene_Id[entrez.row.gene]) != TRUE){
    
  #replace the matched rows with the respective NA spots in the tmp data
  tmp[gene.row,4] <- entrez.row.gene
  
  }
  
    
  #return the data
  return(tmp)
}

#apply the entrez function
new.tmp <- entrez(cna.na.rows,gene.na.rows)
new.tmp[1:5,]
dim(new.tmp)














