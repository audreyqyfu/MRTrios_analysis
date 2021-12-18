# Title: Data preparation and Adding clinical data to find the significant associations 

# Description: Preparing and divide the data by ER test status in order to apply PC and get PCAs  


# Created by Mohamed Megheib

# Date: 05-01-2021

# Last updated: 11-12-2021


#===================================================================================================================================
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

#======================================================methylation================================================================

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

#===========================================Combine expression with meth==========================================================================

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
load("combined.gene.expression.methylation.4.RData")
combined.gene.expression.methylation.4[1:5, 1:6]
dim(combined.gene.expression.methylation.4)

#================================================Adding clinical data===========================================================


clinical.data <- read.delim("brca_tcga_clinical_data.txt", comment.char="#")

#head(clinical.data)
dim(clinical.data)
#names(clinical.data)


ER.clinical.data <- clinical.data[, c(3,28)]

dim(ER.clinical.data)
head(ER.clinical.data)


# changing gene names, i.e TCGA-3C-AAAU-01 to TCGA.3C.AAAU.01

ER.clinical.data$Sample.ID <- gsub("\\-", ".", ER.clinical.data$Sample.ID)
dim(ER.clinical.data)
head(ER.clinical.data)



# matching clinical data with CNA=Meth=gene.express

match.full.clinical <- match(rownames(combined.gene.expression.methylation.4),ER.clinical.data$Sample.ID)
length(match.full.clinical)

match.full.clinical.nona <-match.full.clinical[!is.na(match.full.clinical)]

length(match.full.clinical.nona)


n.ER.clinical.data <- ER.clinical.data[match.full.clinical.nona, ]

dim(n.ER.clinical.data)
rownames(n.ER.clinical.data) <- n.ER.clinical.data$Sample.ID

dim(n.ER.clinical.data)
head(n.ER.clinical.data)



rownames(n.ER.clinical.data)==rownames(combined.gene.expression.methylation.4)

length(rownames(n.ER.clinical.data))==length(rownames(combined.gene.expression.methylation.4))




# Checking dimension 
length(colnames(t.n.ER.clinical.data))==length(colnames(n.n.data_CNA))

n.ER.clinical.data.2 <- n.ER.clinical.data[,-1 ]


combined.gene.expression.methylation.5 <- cbind(data.frame(n.ER.clinical.data.2), combined.gene.expression.methylation.4)

combined.gene.expression.methylation.5[1:4, 1:6]

save(combined.gene.expression.methylation.5, file = "combined.gene.expression.methylation.5.RData")



Negative.combined.gene.expression.methylation.5<- combined.gene.expression.methylation.5[which(combined.gene.expression.methylation.5$n.ER.clinical.data.2== "Negative"), ]
dim(Negative.combined.gene.expression.methylation.5)
Negative.combined.gene.expression.methylation.5[1:4, 1:6]


Positive.combined.gene.expression.methylation.5<- combined.gene.expression.methylation.5[which(combined.gene.expression.methylation.5$n.ER.clinical.data.2== "Positive"), ]
dim(Positive.combined.gene.expression.methylation.5)
Positive.combined.gene.expression.methylation.5[1:4, 1:6]

#==========================================Calculating PCA==================================================================
#Negative
Negative.combined.gene.expression.methylation.pc.2 <- prcomp(Negative.combined.gene.expression.methylation.5[, -1],scale=TRUE)

Negative.summary.combined.gene.expression <- summary(Negative.combined.gene.expression.methylation.pc.2)$importance

write.table(Negative.summary.combined.gene.expression, file = "Negative.summary.combined.gene.expression.txt")

Negative.PCs.combined.gene.expression.methylation.2 <- Negative.combined.gene.expression.methylation.pc.2$x

save(Negative.PCs.combined.gene.expression.methylation.2, file="Negative.PCs.combined.gene.expression.methylation.2.RData")

#Positive
Positive.combined.gene.expression.methylation.pc.2 <- prcomp(Positive.combined.gene.expression.methylation.5[, -1],scale=TRUE)

Positive.summary.combined.gene.expression <- summary(Positive.combined.gene.expression.methylation.pc.2)$importance

write.table(Positive.summary.combined.gene.expression, file = "Positive.summary.combined.gene.expression.txt")

Positive.PCs.combined.gene.expression.methylation.2 <- Positive.combined.gene.expression.methylation.pc.2$x
save(Positive.PCs.combined.gene.expression.methylation.2, file="Positive.PCs.combined.gene.expression.methylation.2.RData")

#===============================================Calculating p-values=======================================================
#Negative 

load("Negative.PCs.combined.gene.expression.methylation.2.RData")
dim(Negative.PCs.combined.gene.expression.methylation.2)
Negative.PCs.combined.gene.expression.methylation.2[1:5, 1:4]

Negative.combined.gene.expression.methylation.5<- combined.gene.expression.methylation.5[which(combined.gene.expression.methylation.5$n.ER.clinical.data.2== "Negative"), ]
dim(Negative.combined.gene.expression.methylation.5)
Negative.combined.gene.expression.methylation.5[1:4, 1:6]

#Use spliting techniques through propagate R Package 
library(propagate)
require(tmvtnorm)
#https://www.rdocumentation.org/packages/propagate/versions/1.0-6/topics/bigcor
cor <- bigcor(Negative.PCs.combined.gene.expression.methylation.2,Negative.combined.gene.expression.methylation.5[, -1] )


Negative.cor.2 <- cor[1: 167, 1: 314707]



Negative.cor.2[1:4, 1:5]

# Calculating t score
#https://courses.lumenlearning.com/introstats1/chapter/testing-the-significance-of-the-correlation-coefficient/
cor.t.score= Negative.cor.2*(sqrt(167-2))/(sqrt(1-Negative.cor.2^2))

cor.t.score[1:4, 1:5]

#https://www.cyclismo.org/tutorial/R/pValues.html
# Calculating p-value

N.cor.t.test.p.value.2=2 * (1 -pt(q = abs(cor.t.score), df = 165))
dim(N.cor.t.test.p.value.2)
N.cor.t.test.p.value.2[1:4, 1:5]

#save(N.cor.t.test.p.value.2, file = "N.cor.t.test.p.value.2.RData")


#Adding key rows to the p-value data sets

N.cor.t.test.p.value.3 <- t(N.cor.t.test.p.value.2)

N.cor.t.test.p.value.3[1:4, 1:5]

#Key rows 
combined.gene.expression.methylation.2[1:3, 1:6]

combined.gene.expression.methylation.2.code <- t(combined.gene.expression.methylation.2[1:3, ])
combined.gene.expression.methylation.2.code[, 1:4]


dim(combined.gene.expression.methylation.2.code)
combined.gene.expression.methylation.2.code[1:4, ]


N.cor.t.test.p.value.3 <- cbind((combined.gene.expression.methylation.2.code),as.data.frame(N.cor.t.test.p.value.3) )
dim(N.cor.t.test.p.value.3)
N.cor.t.test.p.value.3[1:4, 1:5]

colnames(N.cor.t.test.p.value.3) <- row.names(combined.gene.expression.methylation.2)
N.cor.t.test.p.value.3[1:4, 1:5]


save(N.cor.t.test.p.value.3, file = "N.cor.t.test.p.value.3.RData")

#================================================================================

#Positive 

load("Positive.PCs.combined.gene.expression.methylation.2.RData")
dim(Positive.PCs.combined.gene.expression.methylation.2)
Positive.PCs.combined.gene.expression.methylation.2[1:5, 1:4]

Positive.combined.gene.expression.methylation.5<- combined.gene.expression.methylation.5[which(combined.gene.expression.methylation.5$n.ER.clinical.data.2== "Positive"), ]
dim(Positive.combined.gene.expression.methylation.5)
Positive.combined.gene.expression.methylation.5[1:4, 1:6]

#Use spliting techniques through propagate R Package 
library(propagate)
require(tmvtnorm)
#https://www.rdocumentation.org/packages/propagate/versions/1.0-6/topics/bigcor
cor <- bigcor(Positive.PCs.combined.gene.expression.methylation.2,Positive.combined.gene.expression.methylation.5[, -1] )


Positive.cor.2 <- cor[1: 560, 1: 314707]


Positive.cor.2[1:4, 1:5]

# Calculating t score
#https://courses.lumenlearning.com/introstats1/chapter/testing-the-significance-of-the-correlation-coefficient/
cor.t.score= Positive.cor.2*(sqrt(560-2))/(sqrt(1-Positive.cor.2^2))

cor.t.score[1:4, 1:5]

#https://www.cyclismo.org/tutorial/R/pValues.html
# Calculating p-value

P.cor.t.test.p.value.2=2 * (1 -pt(q = abs(cor.t.score), df = 558))
dim(P.cor.t.test.p.value.2)
P.cor.t.test.p.value.2[1:4, 1:5]

#save(N.cor.t.test.p.value.2, file = "N.cor.t.test.p.value.2.RData")


#Adding key rows to the p-value data sets

P.cor.t.test.p.value.3 <- t(P.cor.t.test.p.value.2)

P.cor.t.test.p.value.3[1:4, 1:5]

#Key rows 
combined.gene.expression.methylation.2[1:3, 1:6]

combined.gene.expression.methylation.2.code <- t(combined.gene.expression.methylation.2[1:3, ])
combined.gene.expression.methylation.2.code[, 1:4]


dim(combined.gene.expression.methylation.2.code)
combined.gene.expression.methylation.2.code[1:4, ]


P.cor.t.test.p.value.3 <- cbind((combined.gene.expression.methylation.2.code),as.data.frame(P.cor.t.test.p.value.3) )
dim(P.cor.t.test.p.value.3)
P.cor.t.test.p.value.3[1:4, 1:5]


save(P.cor.t.test.p.value.3, file = "P.cor.t.test.p.value.3.RData")

load("P.cor.t.test.p.value.3.RData")

#=====================================Saving P-values===================================================================
#Negative

load("n.n.expression.4.final.RData")
n.n.expression.2[1:5, 1:4]
dim(n.n.expression.2)

load("N.cor.t.test.p.value.3.RData")
N.cor.t.test.p.value.3[1:5, 1:4]
dim(N.cor.t.test.p.value.3)

Negative.n.n.methylation.3.pvalue <- N.cor.t.test.p.value.3[1:296541, ]
Negative.n.n.expression.2.a <- N.cor.t.test.p.value.3[296542:314707, ]

Negative.n.n.methylation.3.pvalue[1:4, 1:6]
dim(Negative.n.n.methylation.3.pvalue)

Negative.n.n.expression.2.a[1:4, 1:6]
dim(Negative.n.n.expression.2.a)

Negative.n.n.expression.2.a[, 1:3] <- n.n.expression.2[, 1:3]

Negative.n.n.expression.3.pvalue <- Negative.n.n.expression.2.a

save(Negative.n.n.expression.3.pvalue, file = "Negative.n.n.expression.3.pvalue.RData")
save(Negative.n.n.methylation.3.pvalue, file = "Negative.n.n.methylation.3.pvalue.RData")



#Positive

load("n.n.expression.4.final.RData")
n.n.expression.2[1:5, 1:4]
dim(n.n.expression.2)

load("P.cor.t.test.p.value.3.RData")
P.cor.t.test.p.value.3[1:5, 1:4]
dim(P.cor.t.test.p.value.3)

Positive.n.n.methylation.3.pvalue <- P.cor.t.test.p.value.3[1:296541, ]
Positive.n.n.expression.2.a <- P.cor.t.test.p.value.3[296542:314707, ]

Positive.n.n.methylation.3.pvalue[1:4, 1:6]
dim(Positive.n.n.methylation.3.pvalue)

Positive.n.n.expression.2.a[1:4, 1:6]

dim(Positive.n.n.expression.2.a)

Positive.n.n.expression.2.a[, 1:3] <- n.n.expression.2[, 1:3]


Positive.n.n.expression.3.pvalue <- Positive.n.n.expression.2.a

save(Positive.n.n.expression.3.pvalue, file = "Positive.n.n.expression.3.pvalue.RData")
save(Positive.n.n.methylation.3.pvalue, file = "Positive.n.n.methylation.3.pvalue.RData")






