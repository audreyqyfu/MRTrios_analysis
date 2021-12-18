# Title: Data preparation to find the significant associations 

# Description: Preparing and merge the data in order to apply PCA and get PCs  

# Created by Mohamed Megheib

# Date: 05-01-2021

# Last updated: 11-12-2021


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

load("n.n.methylation.4.final.RData")
dim(n.n.methylation.3)

n.n.methylation.3[1:5, 1:4]

# Data containing NAs
missing =n.n.methylation.3[rownames(n.n.methylation.3)[!complete.cases(n.n.methylation.3)], ]

# Data NOT containing NAs
n.n.methylation.3.nona =n.n.methylation.3[rownames(n.n.methylation.3)[complete.cases(n.n.methylation.3)], ]

dim(n.n.methylation.3.nona)

n.n.methylation.5 <- as.data.frame(t(n.n.methylation.3.nona))

n.n.methylation.5[1:5, 1:4]


n.n.methylation.5 <- matrix(c(as.numeric(unlist(n.n.methylation.5))),dim(n.n.methylation.5)[1],dim(n.n.methylation.5)[2])
n.n.methylation.5[1:5, 1:4]
dim(n.n.methylation.5)



#==========================================Gene expression ============================================================

load("n.n.methylation.4.final.RData")
dim(n.n.methylation.3)

load("n.n.expression.4.final.RData")

# Gene expression 
n.n.expression.2[1:5, 1:4]
n.n.expression.3 <- as.data.frame(t(n.n.expression.2))
dim(n.n.expression.3)
n.n.expression.3[1:5, 1:4]

n.n.expression.4 <- n.n.expression.3
n.n.expression.4[1:5, 1:4]
dim(n.n.expression.4)
#unlist the data
n.n.expression.5 <- matrix(c(as.numeric(unlist(n.n.expression.4))),dim(n.n.expression.4)[1],dim(n.n.expression.4)[2])
n.n.expression.5[1:5, 1:4]
dim(n.n.expression.5)

#=====================================================================================================================

#Combine gene expression with methylation data 

combined.gene.expression.methylation.2 <- cbind(n.n.methylation.5,n.n.expression.4)

dim(combined.gene.expression.methylation.2)

combined.gene.expression.methylation.2[1:5, 1:6]

save(combined.gene.expression.methylation.2, file = "combined.gene.expression.methylation.2.RData")

combined.gene.expression.methylation.2[1:5, 1:6]


combined.gene.expression.methylation.3 <- matrix(c(as.numeric(unlist(combined.gene.expression.methylation.2))),dim(combined.gene.expression.methylation.2)[1],dim(combined.gene.expression.methylation.2)[2])
combined.gene.expression.methylation.3[1:5, 1:6]

# Deleting first rows
combined.gene.expression.methylation.4 <- combined.gene.expression.methylation.3[-1:-3, ]
combined.gene.expression.methylation.4[1:5, 1:6]


save(combined.gene.expression.methylation.4, file = "combined.gene.expression.methylation.4.RData")


#============================================================================================================
# Obtaining PCAs
combined.gene.expression.methylation.pc.2 <- prcomp(combined.gene.expression.methylation.4,scale=TRUE)
PCs.combined.gene.expression.methylation.2 <- combined.gene.expression.methylation.pc.2$x


save(PCs.combined.gene.expression.methylation.2, file="PCs.combined.gene.expression.methylation.2.RData")


#=========================================================================================================================
# Reading PCAs
load("PCs.combined.gene.expression.methylaion.2.RData")
dim(PCs.combined.gene.expression.methylation.2)
PCs.combined.gene.expression.methylation.2[1:5, 1:4]

load("combined.gene.expression.methylation.4.RData")
dim(combined.gene.expression.methylation.4)
combined.gene.expression.methylation.4[1:5, 1:4]


#Use spliting techniques through propagate R Package 
library(propagate)
require(tmvtnorm)
#https://www.rdocumentation.org/packages/propagate/versions/1.0-6/topics/bigcor
cor <- bigcor(PCs.combined.gene.expression.methylation.2,combined.gene.expression.methylation.4 )


cor.2 <- cor[1: 770, 1: 314707]


save(cor.2, file = "cor.full.g.m.3.RData")

#write.table(cor, file = "cor.full.g.m.txt")

#Validation

cor.4[1:4, 1:5]

corr.PCs <- corr.test(PCs.combined.gene.expression.methylation[,1:4],data.frame(combined.gene.expression.methylation[, 1:5]),use = 'pairwise.complete.obs')
corr.PCs$r

#================================================================================================

load("cor.full.g.m.3.RData")

cor.2[1:4, 1:5]

# Calculating t score
#https://courses.lumenlearning.com/introstats1/chapter/testing-the-significance-of-the-correlation-coefficient/
cor.t.score= cor.2*(sqrt(770-2))/(sqrt(1-cor.2^2))

cor.t.score[1:4, 1:5]

#https://www.cyclismo.org/tutorial/R/pValues.html
# Calculating p-value

cor.t.test.p.value.2=2 * (1 -pt(q = abs(cor.t.score), df = 768))
dim(cor.t.test.p.value.2)
cor.t.test.p.value.2[1:4, 1:5]

save(cor.t.test.p.value.2, file = "cor.t.test.p.value.2.RData")


#Adding key rows to the p-value data sets

cor.t.test.p.value.3 <- t(cor.t.test.p.value.2)

cor.t.test.p.value.3[1:4, 1:5]

#Key rows 
combined.gene.expression.methylation.2[1:3, 1:6]

combined.gene.expression.methylation.2.code <- t(combined.gene.expression.methylation.2[1:3, ])
combined.gene.expression.methylation.2.code[, 1:4]


dim(combined.gene.expression.methylation.2.code)
combined.gene.expression.methylation.2.code[1:4, ]


cor.t.test.p.value.3 <- cbind((combined.gene.expression.methylation.2.code),as.data.frame(cor.t.test.p.value.3) )
dim(cor.t.test.p.value.3)
cor.t.test.p.value.3[1:4, 1:5]

colnames(cor.t.test.p.value.3) <- row.names(combined.gene.expression.methylation.2)


save(cor.t.test.p.value.3, file = "cor.t.test.p.value.3.RData")


#======================================================================================================
# Calculating q-values 


load("cor.t.test.p.value.3.RData")
dim(cor.t.test.p.value.3)
cor.t.test.p.value.3[1:4, 1:5]


cor.t.test.p.value.t <- t(cor.t.test.p.value.3)
dim(cor.t.test.p.value.t)
cor.t.test.p.value.t[1:4, 1:5]

All.Pvalues <- list()
List.significant.gene.expression.methylation.TCGA.5 <- list()

All.qvalues <- matrix(0,dim(cor.t.test.p.value)[1], dim(cor.t.test.p.value)[2])

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

sum(All.qvalues.Significant[770, ])

save(List.significant.gene.expression.methylation.TCGA.5, file="List.significant.gene.expression.methylation.TCGA.5.RData")


#========================================================================================================
#

load("cor.t.test.p.value.3.RData")

n.n.methylation.3.pvalue <- cor.t.test.p.value.3[1:296541, ]
n.n.expression.2.a <- cor.t.test.p.value.3[296542:314707, ]

n.n.methylation.3.a[1:4, 1:6]
dim(n.n.methylation.3.a)

n.n.expression.2.a[1:4, 1:6]
n.n.expression.2[1:4, 1:6]

dim(n.n.methylation.3)

Negative.n.n.expression.2.a[, 1:3] <- n.n.expression.2[, 1:3]

data_apply <- n.n.expression.2.a[, -1:-3]

# Replicate original data
data_new <- n.n.expression.2        

# Replace specific columns
data_new[ , colnames(data_new) %in% colnames(data_apply)] <- data_apply  
data_new[1:4, 1:5] 

n.n.expression.2.pvalue <- data_new

save(n.n.expression.2.pvalue, file = "n.n.expression.2.pvalue.RData")
save(n.n.methylation.3.pvalue, file = "n.n.methylation.3.pvalue.RData")


#===========================================================================================================

