# Title: Second independent approach to identify significant associations 

# Description: the p values of gene expression+ methylation are combined then apply qvalue methods for each individual pc.  

# Created by Mohamed Megheib

# Date: 11-01-2021

# Last updated: 12-12-2021


#======================================================================================================
# Reading methylation  and gene expression 
load("n.n.methylation.4.final.RData")
load("n.n.expression.4.final.RData")

# Reading P-values 

load("cor.t.test.p.value.3.RData")
dim(cor.t.test.p.value.3)
cor.t.test.p.value.3[1:4, 1:5]

#Extract identifiers 
code.cor.t.test.p.value.3 <- cor.t.test.p.value.3[, 1:3]

cor.t.test.p.value.t <- t(cor.t.test.p.value)
dim(cor.t.test.p.value.t)
cor.t.test.p.value.t[1:4, 1:5]

All.Pvalues <- list()
List.significant.gene.expression.methylation.TCGA.5 <- list()

All.qvalues.Significant <- matrix(0,dim(cor.t.test.p.value.3)[1], dim(cor.t.test.p.value.3)[2]-3)

# calculate associated PCs 

cor.t.test.p.value.3 <- cor.t.test.p.value.3[ ,-1:-3]

for (m in 1:dim(cor.t.test.p.value.3)) {
  # The p values
  
  Pvalues <- cor.t.test.p.value.3[, m]
  Pvalues.nona <- Pvalues[!is.na(Pvalues)]
  qobj <- qvalue(Pvalues.nona, fdr.level=0.10) 
  
  # Significant associations
  All.qvalues.Significant[,m] <- qobj$significant
  ## List.significant.gene.expression.methylation.TCGA.5[[m]] <- which(Significant.asso,useNames = TRUE)
  
  print(m)
}


All.qvalues.Significant[1:5, 1:6]

# Adding identifiers 
All.qvalues.Significant.2 <- cbind(code.cor.t.test.p.value.3,All.qvalues.Significant )
All.qvalues.Significant.2[1:5, 1:6]
dim(All.qvalues.Significant.2)



#========================================================================================================
# Separate the gene expression from methylation

n.n.methylation.4.pvalue <- All.qvalues.Significant.2[1:296541, ]
n.n.expression.4.pvalue <- All.qvalues.Significant.2[296542:314707, ]

n.n.methylation.3[1:4, 1:6]
n.n.methylation.4.pvalue[1:4, 1:6]
n.n.expression.2[1:4, 1:6]
n.n.expression.4.pvalue[1:4, 1:6]

n.n.expression.4.pvalue[, 1:3] <- n.n.expression.2[, 1:3]


save(n.n.expression.4.pvalue, file = "n.n.expression.4.pvalue.RData")
save(n.n.methylation.4.pvalue, file = "n.n.methylation.4.pvalue.RData")

#======================================================================================

#Constructing trios 

n.n.methylation.3 <-  n.n.methylation.4.pvalue
n.n.methylation.3[1:5, 1:6]
dim(n.n.methylation.3)
n.n.expression.2 <-  n.n.expression.4.pvalue
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
    
    #print(exp.meth[1:6, ])
    
    dim(exp.meth)
    exp.meth.1 <- exp.meth[-1:-3,]
    
    exp.meth.1[1:4, ]
    
    exp.meth.2 <- matrix(as.numeric(unlist(exp.meth.1)), dim(exp.meth.1)[1], 2)
    
    exp.meth.2[1:4, ]
    
    dim(exp.meth.2)
    
    
    exp.meth.3 <- exp.meth.2[, 1]+exp.meth.2[, 2]
    
    exp.meth.4 <- sum(exp.meth.3 !=0)
    
    write.table(t(c(unlist(exp.meth[1, 2]), unlist(exp.meth[2, 2]),unlist(exp.meth[3, 2]), exp.meth.4 )), file="Indep.setup.TCGA.8.2nd.txt", sep = "\t", row.names = FALSE,col.names = FALSE , append =TRUE,quote=FALSE)
    
    
  }
  print(length(full.range)-k)
  print("from 18890")
}



#===========================================Negative===============================================================================

# Reading P-values 
load("N.cor.t.test.p.value.3.RData")
N.cor.t.test.p.value.3[1:5, 1:4]
dim(N.cor.t.test.p.value.3)
code.N.cor.t.test.p.value.3 <- N.cor.t.test.p.value.3[, 1:3]

N.cor.t.test.p.value.3[1:4, 1:5]

All.qvalues.Significant <- matrix(0,dim(N.cor.t.test.p.value.3)[1], dim(N.cor.t.test.p.value.3)[2]-3)

dim(All.qvalues.Significant)

# calculate associated PCs 

cor.t.test.p.value.3 <- N.cor.t.test.p.value.3[ ,-1:-3]

for (m in 1:dim(cor.t.test.p.value.3)) {
  # The p values
  
  Pvalues <- cor.t.test.p.value.3[, m]
  Pvalues.nona <- Pvalues[!is.na(Pvalues)]
  qobj <- qvalue(Pvalues.nona, fdr.level=0.10) 
  
  # Significant associations
  All.qvalues.Significant[,m] <- qobj$significant
  ## List.significant.gene.expression.methylation.TCGA.5[[m]] <- which(Significant.asso,useNames = TRUE)
  
  print(m)
}


All.qvalues.Significant[1:5, 1:6]


All.qvalues.Significant.2 <- cbind(code.N.cor.t.test.p.value.3,All.qvalues.Significant )
All.qvalues.Significant.2[1:5, 1:6]
dim(All.qvalues.Significant.2)


load("n.n.expression.4.final.RData")

# Gene expression 
n.n.expression.2[1:5, 1:4]

N.n.n.methylation.4.qvalue <- All.qvalues.Significant.2[1:296541, ]
N.n.n.expression.4.qvalue <- All.qvalues.Significant.2[296542:314707, ]


N.n.n.methylation.4.qvalue[1:4, 1:6]


n.n.expression.2[1:4, 1:6]
N.n.n.expression.4.qvalue[1:4, 1:6]


N.n.n.expression.4.qvalue[, 1:3] <- n.n.expression.2[, 1:3]

#===================================================================

n.n.methylation.3 <-  N.n.n.methylation.4.qvalue
n.n.methylation.3[1:5, 1:6]
dim(n.n.methylation.3)
n.n.expression.2 <-  N.n.n.expression.4.qvalue
n.n.expression.2[1:5, 1:6]
dim(n.n.expression.2)

full.range <- unique(n.n.methylation.3$Entrez_Gene_Id)



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
    
    #print(exp.meth[1:6, ])
    
    dim(exp.meth)
    exp.meth.1 <- exp.meth[-1:-3,]
    
    exp.meth.1[1:4, ]
    
    exp.meth.2 <- matrix(as.numeric(unlist(exp.meth.1)), dim(exp.meth.1)[1], 2)
    
    exp.meth.2[1:4, ]
    
    dim(exp.meth.2)
    
    
    exp.meth.3 <- exp.meth.2[, 1]+exp.meth.2[, 2]
    
    exp.meth.4 <- sum(exp.meth.3 !=0)
    
    write.table(t(c(unlist(exp.meth[1, 2]), unlist(exp.meth[2, 2]),unlist(exp.meth[3, 2]), exp.meth.4 )), file="Indep.setup.TCGA.8.2nd.Negative.txt", sep = "\t", row.names = FALSE,col.names = FALSE , append =TRUE,quote=FALSE)
    
    
  }
  print(length(full.range)-k)
  print("from 18890")
}



#================================================Positive=================================================

# Reading P-values 

load("P.cor.t.test.p.value.3.RData")
P.cor.t.test.p.value.3[1:5, 1:4]
dim(P.cor.t.test.p.value.3)

code.P.cor.t.test.p.value.3 <- P.cor.t.test.p.value.3[, 1:3]

P.cor.t.test.p.value.3[1:4, 1:5]

All.qvalues.Significant <- matrix(0,dim(P.cor.t.test.p.value.3)[1], dim(P.cor.t.test.p.value.3)[2]-3)

dim(All.qvalues.Significant)

# calculate associated PCs 

cor.t.test.p.value.3 <- P.cor.t.test.p.value.3[ ,-1:-3]

for (m in 1:dim(cor.t.test.p.value.3)) {
  # The p values
  
  Pvalues <- cor.t.test.p.value.3[, m]
  Pvalues.nona <- Pvalues[!is.na(Pvalues)]
  qobj <- qvalue(Pvalues.nona, fdr.level=0.10) 
  
  # Significant associations
  All.qvalues.Significant[,m] <- qobj$significant
  ## List.significant.gene.expression.methylation.TCGA.5[[m]] <- which(Significant.asso,useNames = TRUE)
  
  print(m)
}


All.qvalues.Significant[1:5, 1:6]


All.qvalues.Significant.2 <- cbind(code.P.cor.t.test.p.value.3,All.qvalues.Significant )
All.qvalues.Significant.2[1:5, 1:6]
dim(All.qvalues.Significant.2)


load("n.n.expression.4.final.RData")

# Gene expression 
n.n.expression.2[1:5, 1:4]

P.n.n.methylation.4.qvalue <- All.qvalues.Significant.2[1:296541, ]
P.n.n.expression.4.qvalue <- All.qvalues.Significant.2[296542:314707, ]


P.n.n.methylation.4.qvalue[1:4, 1:6]


n.n.expression.2[1:4, 1:6]
P.n.n.expression.4.qvalue[1:4, 1:6]


P.n.n.expression.4.qvalue[, 1:3] <- n.n.expression.2[, 1:3]

#===================================================================

n.n.methylation.3 <-  P.n.n.methylation.4.qvalue
n.n.methylation.3[1:5, 1:6]
dim(n.n.methylation.3)
n.n.expression.2 <-  P.n.n.expression.4.qvalue
n.n.expression.2[1:5, 1:6]
dim(n.n.expression.2)

full.range <- unique(n.n.methylation.3$Entrez_Gene_Id)



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
    
    #print(exp.meth[1:6, ])
    
    dim(exp.meth)
    exp.meth.1 <- exp.meth[-1:-3,]
    
    exp.meth.1[1:4, ]
    
    exp.meth.2 <- matrix(as.numeric(unlist(exp.meth.1)), dim(exp.meth.1)[1], 2)
    
    exp.meth.2[1:4, ]
    
    dim(exp.meth.2)
    
    
    exp.meth.3 <- exp.meth.2[, 1]+exp.meth.2[, 2]
    
    exp.meth.4 <- sum(exp.meth.3 !=0)
    
    write.table(t(c(unlist(exp.meth[1, 2]), unlist(exp.meth[2, 2]),unlist(exp.meth[3, 2]), exp.meth.4 )), file="Indep.setup.TCGA.8.2nd.Positive.txt", sep = "\t", row.names = FALSE,col.names = FALSE , append =TRUE,quote=FALSE)
    
    
  }
  print(length(full.range)-k)
  print("from 18890")
}






















