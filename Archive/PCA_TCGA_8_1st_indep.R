# Title: First independent approach to identify significant associations 

# Description: Constructing trio from gene expresssion and methylation p values, then apply qvalue methods for each trio.  

# Created by Mohamed Megheib

# Date: 10-01-2021

# Last updated: 12-12-2021


#======================================================================================================
# Library
library(stringi,lib="/mnt/ceph/megheib/Rpackages") #need for stri_sub ()
library(stringr,lib="/mnt/ceph/megheib/Rpackages")
library(MRPC,lib="/mnt/ceph/megheib/Rpackages")
library(fastcluster,lib="/mnt/ceph/megheib/Rpackages")
library(WGCNA,lib="/mnt/ceph/megheib/Rpackages")
library(psych,lib="/mnt/ceph/megheib/Rpackages")
library(MatrixEQTL,lib="/mnt/ceph/megheib/Rpackages")

library(MRPC,lib="/mnt/ceph/megheib/Rpackages")

library(mice)

require(dplyr)

#========================================================================================================

#Loading p-values of meth and gene expression. 

load("n.n.expression.3.pvalue.RData")
load("n.n.methylation.3.pvalue.RData")

n.n.methylation.3 <-  n.n.expression.3.pvalue
n.n.methylation.3[1:5, 1:6]
dim(n.n.methylation.3)
n.n.expression.2 <- n.n.expression.3.pvalue
n.n.expression.2[1:5, 1:6]
dim(n.n.expression.2)

full.range <- unique(n.n.methylation.3$Entrez_Gene_Id)


List.significant.gene.expression.methylation.TCGA.5 <- list()

for (k in 1:length(full.range)) {

  exp.meth <- matrix(0)
  exp <- NULL
  meth <- NULL
  
  exp[k] <- paste0(full.range[k])
  meth[k] <- paste0(full.range[k])
  
  # Adding filters to skip when the data do not match
  
  if(dim(n.n.expression.2[which(n.n.expression.2$Entrez_Gene_Id==exp[k]),])[1]==0) next
  
  # defining matrices to save results 
  
  # Creating a loop to include repeated cases 
  
  for (j in 1: dim(n.n.methylation.3[which(n.n.methylation.3$Entrez_Gene_Id==meth[k]),])[1]) {

 
    if(dim(n.n.methylation.3[which(n.n.methylation.3$Entrez_Gene_Id==meth[k]),][j, ])[1]==0) next
    
    exp.meth<-matrix(c(n.n.expression.2[which(n.n.expression.2$Entrez_Gene_Id==exp[k]),],
                           n.n.methylation.3[which(n.n.methylation.3$Entrez_Gene_Id==meth[k]),][j, ]), dim(n.n.expression.2)[2] , 2)
    
   # print(exp.meth[1:6, ])
    
    dim(exp.meth)
    exp.meth.1 <- exp.meth[-1:-3,]
    
    exp.meth.1[1:4, ]  
    

    exp.meth.2 <- matrix(as.numeric(unlist(exp.meth.1)), dim(exp.meth.1)[1], 2)
    
    exp.meth.2[1:4, ]
     
    dim(exp.meth.2)
  
    
    # calculate associated PCs 
    
      # The p values
      Pvalues <- exp.meth.2
      Pvalues.nona <- Pvalues[!is.na(Pvalues)]
      #All.Pvalues [[m]] <- Pvalues.nona
      qobj <- qvalue(Pvalues.nona, fdr.level=0.10) 
      
      # Significant associations
      Significant.asso <- qobj$significant
      
      
      List.significant.gene.expression.methylation.TCGA.5[[j]] <- which(Significant.asso,useNames = TRUE)
         
      write.table(t(c(unlist(exp.meth[1, 2]), unlist(exp.meth[2, 2]),unlist(exp.meth[3, 2]),length(unique(List.significant.gene.expression.methylation.TCGA.5[[j]])) , sort(unique(List.significant.gene.expression.methylation.TCGA.5[[j]])) )), file="Indep.setup.TCGA.8.sort.Positive.txt", sep = "\t", row.names = FALSE,col.names = FALSE , append =TRUE,quote=FALSE)
      
      
    
  }
 # print(length(full.range)-k)
  #print("from 18890")
  }

for(i in 1:length(List.significant.gene.expression.methylation.TCGA.5[[j]])){
  if(List.significant.gene.expression.methylation.TCGA.5[[j]][i]>170) {
    List.significant.gene.expression.methylation.TCGA.5[[j]][i]= List.significant.gene.expression.methylation.TCGA.5[[j]][i]-170}
  else{
    List.significant.gene.expression.methylation.TCGA.5[[j]][i]=List.significant.gene.expression.methylation.TCGA.5[[j]][i]
  }
}







