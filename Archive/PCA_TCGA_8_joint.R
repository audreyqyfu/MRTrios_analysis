# Title: Joint approach to identify significant associations 

# Description: Constructing trio from gene expresssion and methylation data, then use p-values of linear regression, 

#then apply qvalue methods for each trio.  

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

#==========================================methylation================================================================

# Reading methylation data 

load("n.n.methylation.4.final.RData")
dim(n.n.methylation.3)

n.n.methylation.3[1:5, 1:4]

# Data NOT containing NAs
n.n.methylation.3.nona =n.n.methylation.3[rownames(n.n.methylation.3)[complete.cases(n.n.methylation.3)], ]

dim(n.n.methylation.3.nona)

n.n.methylation.3 <- n.n.methylation.3.nona

n.n.methylation.3[1:5, 1:4]

#==========================================Gene expression ============================================================
# Reading gene expression data 

load("n.n.expression.4.final.RData")

# Gene expression 
n.n.expression.2[1:5, 1:4]

#============================================================================================================
# pca matrix

load("PCs.combined.gene.expression.methylation.2.RData")
dim(PCs.combined.gene.expression.methylation.2)
PCs.combined.gene.expression.methylation.2[1:5, 1:4]


full.range <- unique(n.n.methylation.3$Entrez_Gene_Id)
length(full.range)
full.range [1:5]


#n.n.expression.2[which(n.n.expression.2$Entrez_Gene_Id=="10357"),1:4]


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
    
    print(exp.meth[1:6, ])
    
    
    exp.meth.1 <- exp.meth[-1:-3,]
    
    exp.meth.1[1:4, ]  
    
    exp.meth.2 <- matrix(as.numeric(unlist(exp.meth.1)), dim(exp.meth.1)[1], 2)
    
    exp.meth.2[1:4, ]
    
    dim(exp.meth.2)
    
    
    Pvalues <- matrix(0)
    
    for (i in 1: dim(PCs.combined.gene.expression.methylation.2)[1]) {
      
    l <- lm(PCs.combined.gene.expression.methylation.2[, i]~ exp.meth.2[, 1]+exp.meth.2[, 2])
    
    f = summary(l)$fstatistic
    
    P= pf(f[1],f[2],f[3],lower.tail=F)
    
    Pvalues[i] <- P
    
    
    }
    
    #All.Pvalues <- list()
    # calculate associated PCs 
    # The p values
    Pvalues.nona <- Pvalues[!is.na(Pvalues)]
    #All.Pvalues [[m]] <- Pvalues.nona
    qobj <- qvalue(Pvalues.nona, fdr.level=0.10) 
    
    # Significant associations
    Significant.asso <- qobj$significant
    List.significant.gene.expression.methylation.TCGA.5[[j]] <- which(Significant.asso,useNames = TRUE)
    

    write.table(t(c(unlist(exp.meth[1, 2]), unlist(exp.meth[2, 2]),unlist(exp.meth[3, 2]),length(List.significant.gene.expression.methylation.TCGA.5[[j]]) , List.significant.gene.expression.methylation.TCGA.5[[j]] )), file="Joint.setup.TCGA.6.txt", sep = "\t", row.names = FALSE,col.names = FALSE , append =TRUE,quote=FALSE)

    
    print(length(full.range)-k)
    print("from 18890")
    
    
     
  }
}



