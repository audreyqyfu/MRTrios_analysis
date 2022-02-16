
#Title: TCGA data Preparation

# Description: cleaning the data, adding and updating entrez ids, matching, dividing the data by ER status.

# Created by Mohamed Megheib

# Date: 05-01-2021

# Last updated: 11-12-2021


#===================================================================================================================================

library(MRPC,lib="/mnt/ceph/megheib/Rpackages")

require(dplyr)


#===================================================================================================================================

# Reading methylation data 

load("MM.t.6.RData")

MM.t.6[1:5, 1:6]

dim(MM.t.6)

#Split delimited strings in a column and insert as new rows

MM.t.6 <- separate_rows(MM.t.6,Hugo_Symbol, sep = ";")

MM.t.6[1:4, 1:5]

dim(MM.t.6)

table(MM.t.6$Hugo_Symbol)[1:4]


# delete rows with no names 
MM.t.7 <- MM.t.6[!MM.t.6$Hugo_Symbol=="",]
dim(MM.t.7)
MM.t.7[1:5, 1:6]



#====================================================================================================================================


# Reading entrez ID Lsit

ntrez_id=read.csv("Genes_entrezID_h38_ensembl_biomart.txt")


# Using package "org.Hs.eg.db" to complete entrez id

library(org.Hs.eg.db)

hs <- org.Hs.eg.db

my.symbols <- MM.t.7$Hugo_Symbol

my.symbols[1:15]

ntrez_id <- select(hs, 
               keys = my.symbols,
               columns = c("ENTREZID", "SYMBOL"),
               keytype = "SYMBOL")

dim(ntrez_id)
colnames(ntrez_id) <- c("Hugo_Symbol", "Entrez_Gene_Id")

ntrez_id[1:5, ]


# insert entrez id inside methylation 


full.data <- MM.t.7$Hugo_Symbol

full.data[1:6]

id <- length(full.data)

for (i in 1:id) {
  
  meth <- NULL
  
  meth[i] <- paste0(full.data[i])
  
 # print( meth[i])
  
 MM.t.7$Entrez_Gene_Id[i] <- matrix(c(((ntrez_id[which(ntrez_id$Hugo_Symbol==meth[i]),]))), 1, 2)[, 2]
  #print(MM.t.7$Entrez_Gene_Id[i])
  #print(matrix(c(((ntrez_id[which(ntrez_id$Hugo_Symbol==meth[i]),]))), 1, 2)[, 2])
}


ntrez_id[which(ntrez_id$Hugo_Symbol=="ACTN1"),]

matrix(c(((ntrez_id[which(ntrez_id$Hugo_Symbol=="ACTN1"),]))), 1, 2)[, 2]

MM.t.7[1:4 , c(1:3, 899)]

MM.t.7 <- MM.t.7 %>% relocate(Entrez_Gene_Id, .before = Genom)

MM.t.7[1:5, 1:6]

#save(MM.t.7, file = "MM.t.7.7.RData")


MM.t.7[which(MM.t.7$Entrez_Gene_Id==5577), 1:5]


#========================================================================================================================
#CNA
data_CNA=read.delim("data_CNA.txt")
head(data_CNA)[1:3]
dim(data_CNA)

data_CNA[which(data_CNA$Entrez_Gene_Id=="647042"), 1:5]

#Extracting a list of probes with entrez id 

prob.geneid.list <- MM.t.7[, c(4, 6)]

dim(prob.geneid.list)

prob.geneid.list[1:5,]

# Insert probe id inside CNA 

full.data <- data_CNA$Entrez_Gene_Id

id <- length(full.data)

for (i in 1:id) {
  
  meth <- NULL
  
  meth[i] <- paste0(full.data[i])
  
  data_CNA$Probe.ID[i] <- matrix(c(unlist((prob.geneid.list[which(prob.geneid.list$Probe.ID==meth[i]),]))), 1, 2)[, 1]
  
}


data_CNA[1:4 , c(1:3, 899)]

CNA.t.7 <- data_CNA %>% relocate(Probe.ID, .before = Hugo_Symbol)

CNA.t.7[1:5, 1:6]

#==========================================================================================================================================

#expression
expression <- read.delim("data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt")
head(expression)[1:5]

dim(expression)


# insert probe id inside expression 

full.data <- expression$Entrez_Gene_Id

id <- length(full.data)

for (i in 1:id) {
  
  meth <- NULL
  
  meth[i] <- paste0(full.data[i])
  
  expression$Probe.ID[i] <- matrix(c(unlist((prob.geneid.list[which(prob.geneid.list$Probe.ID==meth[i]),]))), 1, 2)[, 1]
  
}


expression[1:4 , c(1:3, 899)]

EM.t.7 <- expression %>% relocate(Probe.ID, .before = Hugo_Symbol)

EM.t.7[1:5, 1:6]


dim(CNA.t.7)
dim(EM.t.7)
dim(MM.t.7)

CNA.t.7[1:3, 1:4]
EM.t.7[1:3, 1:4]
MM.t.7[1:3, 1:6]

table(CNA.t.7$Hugo_Symbol)[1:4]
table(EM.t.7$Hugo_Symbol)[1:4]
table(MM.t.7$Hugo_Symbol)[1:4]

#================================================================================================================================

# Matching individual together in Meth, gene expression and CNA

align.col=match(colnames(MM.t.7),colnames(CNA.t.7))

MM.t.8=MM.t.7[, -which(is.na(align.col))]

dim(MM.t.8)

MM.t.8[1:5,1:5]


CNA.t.7=CNA.t.7[, na.omit(align.col)]

dim(CNA.t.7)
CNA.t.7[1:5,1:5]


align.col.2=match(colnames(CNA.t.7),colnames(EM.t.7))


EM.t.10=EM.t.7[, na.omit(align.col.2)]
dim(EM.t.10)
EM.t.10[1:5,1:5]


align.col.3=match(colnames(EM.t.10),colnames(CNA.t.7))

CNA.t.10=CNA.t.7[, na.omit(align.col.3)]

dim(CNA.t.10)
CNA.t.10[1:5,1:5]


align.col.4=match(colnames(CNA.t.10),colnames(MM.t.8))

MM.t.10=MM.t.8[, na.omit(align.col.4)]

dim(MM.t.10)
MM.t.10[1:5,1:5]

#===================================================================================================================================

CNA.t.10.unique<- CNA.t.10[!duplicated(CNA.t.10[, 3]), ]
#EM.t.10.unique <- EM.t.10[!duplicated(EM.t.10[, 3]), ]

dim(CNA.t.10.unique)

CNA.t.10.unique[1:5, 1:6]

save(CNA.t.10.unique, file ="CNA.t.11.unique.RData")
save(EM.t.10, file ="EM.t.11.RData")
save(MM.t.10, file ="MM.t.11.RData")


#====================================================================================================================================

save(CNA.t.10,file="/mnt/ceph/megheib/Cancer_Genomics/CNA.t.10.RData")
save(EM.t.10,file="/mnt/ceph/megheib/Cancer_Genomics/EM.t.10.RData")
save(MM.t.10,file="/mnt/ceph/megheib/Cancer_Genomics/MM.t.10.RData")

#===================================================================================================================================

load("CNA.t.11.unique.RData")
load("EM.t.11.RData")
load("MM.t.11.RData")

dim(CNA.t.10.unique)
dim(EM.t.10)
dim(MM.t.10)


CNA.t.10.unique[1:5, 1:5]
EM.t.10[1:5, 1:5]
MM.t.10[1:5, 1:5]


MM.t.10[which(MM.t.10$Hugo_Symbol=="RBL2"), 1:5]
EM.t.10[which(EM.t.10$Hugo_Symbol=="RBL2"), 1:5]
CNA.t.10[which(CNA.t.10$Hugo_Symbol=="RBL2"), 1:5]



MM.t.10[which(MM.t.10$Entrez_Gene_Id==5934), 1:5]
EM.t.10[which(EM.t.10$Entrez_Gene_Id==5934), 1:5]
CNA.t.10[which(CNA.t.10$Entrez_Gene_Id==5934), 1:5]


CNA.t.10.unique[which(CNA.t.10.unique$Entrez_Gene_Id==400798), 1:5]
n.n.data_CNA.2[which(n.n.data_CNA.2$Entrez_Gene_Id==645425), 1:5]
n.n.data_CNA.2[which(n.n.data_CNA.2$Entrez_Gene_Id==728695), 1:5]



CNA.t.10.unique<- CNA.t.10[!duplicated(CNA.t.10[, 3]), ]

dim(CNA.t.10.unique)

save(CNA.t.10.unique, file = "CNA.t.10.unique.RData")

CNA.t.10.unique[1:5, 1:6]


save(CNA.t.10.unique,file="/mnt/ceph/megheib/Cancer_Genomics/CNA.t.10.unique.RData")
save(EM.t.10,file="/mnt/ceph/megheib/Cancer_Genomics/EM.t.10.RData")
save(MM.t.10,file="/mnt/ceph/megheib/Cancer_Genomics/MM.t.10.RData")




#==============================================================================================================================

# Adding ER_status

#==============================================================================================================================
# Reading CNA and Gene expression 

load("CNA.t.9.unique.RData")
load("EM.t.9.unique.RData")
load("MM.t.9.up.RData")


n.n.data_CNA <- CNA.t.9.unique
n.n.expression <- EM.t.9.unique
n.n.methylation <- MM.t.9
n.n.methylation.na <- n.n.methylation.2[is.na(n.n.methylation.2$Entrez_Gene_Id), ]


dim(n.n.data_CNA)
dim(n.n.methylation)
dim(n.n.expression)

n.n.data_CNA[1:4, 1:5]
n.n.expression[1:4, 1:5]
n.n.methylation[1:4, 1:5]



#===================================================================================================================================

######################Adding demographics to the main data##########


clinical.data <- read.delim("brca_tcga_clinical_data.txt", comment.char="#")

#head(clinical.data)
dim(clinical.data)
#names(clinical.data)


ER.clinical.data <- clinical.data[, c(3,28)]

#dim(ER.clinical.data)
head(ER.clinical.data)


# changing gene names, i.e TCGA-3C-AAAU-01 to TCGA.3C.AAAU.01

ER.clinical.data$Sample.ID <- gsub("\\-", ".", ER.clinical.data$Sample.ID)

#head(ER.clinical.data)


# matching clinical data with CNA=Meth=gene.express

match.full.clinical <- match(colnames(n.n.data_CNA),ER.clinical.data$Sample.ID)
length(match.full.clinical)

match.full.clinical.nona <-match.full.clinical[!is.na(match.full.clinical)]

length(match.full.clinical.nona)


n.ER.clinical.data <- ER.clinical.data[match.full.clinical.nona, ]

dim(n.ER.clinical.data)

# Rotating clinical data to match with others 
t.n.ER.clinical.data <- t(n.ER.clinical.data)

t.n.ER.clinical.data <- as.data.frame(t.n.ER.clinical.data)

colnames(t.n.ER.clinical.data) <- n.ER.clinical.data$Sample.ID


t.n.ER.clinical.data <- t.n.ER.clinical.data[-1,]

dim(t.n.ER.clinical.data)

#t.n.ER.clinical.data[, 1:5]


t.n.ER.clinical.data$Probe.ID <- "ER"
t.n.ER.clinical.data$Hugo_Symbol <- " "
t.n.ER.clinical.data$Entrez_Gene_Id <- " "


t.n.ER.clinical.data <- t.n.ER.clinical.data %>% relocate(Entrez_Gene_Id, .before = TCGA.A2.A0YK.01)

t.n.ER.clinical.data <- t.n.ER.clinical.data %>% relocate(Hugo_Symbol, .before = Entrez_Gene_Id)

t.n.ER.clinical.data <- t.n.ER.clinical.data %>% relocate(Probe.ID, .before = Hugo_Symbol)

t.n.ER.clinical.data[, 1:5]

dim(t.n.ER.clinical.data)

# Checking dimension 
length(colnames(t.n.ER.clinical.data))==length(colnames(n.n.data_CNA))

#===================================================================================================================================
#CNA

# Combine clinical data with main CNA
n.n.n.data_CNA <- rbind(n.n.data_CNA, t.n.ER.clinical.data)
dim(n.n.n.data_CNA)
n.n.n.data_CNA[1:5, 1:6]

# Identify ER column 
which(n.n.n.data_CNA$Probe.ID=="ER")

# Choose negative 

Negative.n.n.n.data_CNA<- n.n.n.data_CNA[, which(n.n.n.data_CNA[which(n.n.n.data_CNA$Probe.ID=="ER"), ]== "Negative")]
Negative.n.n.n.data_CNA$Probe.ID <- n.n.n.data_CNA$Probe.ID
Negative.n.n.n.data_CNA$Hugo_Symbol <- n.n.n.data_CNA$Hugo_Symbol
Negative.n.n.n.data_CNA$Entrez_Gene_Id <- n.n.n.data_CNA$Entrez_Gene_Id


Negative.n.n.n.data_CNA <- Negative.n.n.n.data_CNA %>% relocate(Entrez_Gene_Id, .before = TCGA.A2.A3XT.01)
Negative.n.n.n.data_CNA <- Negative.n.n.n.data_CNA %>% relocate(Hugo_Symbol, .before = Entrez_Gene_Id)
Negative.n.n.n.data_CNA <- Negative.n.n.n.data_CNA %>% relocate(Probe.ID, .before = Hugo_Symbol)

dim(Negative.n.n.n.data_CNA)

Negative.n.n.n.data_CNA[1:4, 1:5]


#G.expression
# Combine clinical data with main CNA
n.n.n.expression <- rbind(n.n.expression, t.n.ER.clinical.data)
dim(n.n.n.expression)
n.n.n.expression[1:5, 1:6]

# Identify ER column 
which(n.n.n.expression$Probe.ID=="ER")

# Choose negative 

Negative.n.n.n.expression<- n.n.n.expression[, which(n.n.n.expression[which(n.n.n.expression$Probe.ID=="ER"), ]== "Negative")]
Negative.n.n.n.expression$Probe.ID <- n.n.n.expression$Probe.ID
Negative.n.n.n.expression$Hugo_Symbol <- n.n.n.expression$Hugo_Symbol
Negative.n.n.n.expression$Entrez_Gene_Id <- n.n.n.expression$Entrez_Gene_Id


Negative.n.n.n.expression <- Negative.n.n.n.expression %>% relocate(Entrez_Gene_Id, .before = TCGA.A2.A3XT.01)
Negative.n.n.n.expression <- Negative.n.n.n.expression %>% relocate(Hugo_Symbol, .before = Entrez_Gene_Id)
Negative.n.n.n.expression <- Negative.n.n.n.expression %>% relocate(Probe.ID, .before = Hugo_Symbol)

dim(Negative.n.n.n.expression)

Negative.n.n.n.expression[1:4, 1:5]


#methylation
#Combine clinical data with main CNA
n.n.n.methylation<- rbind(n.n.methylation, t.n.ER.clinical.data)
dim(n.n.n.methylation)
n.n.n.methylation[1:5, 1:6]

# Identify ER column 
which(n.n.n.methylation$Probe.ID=="ER")

# Choose negative 

Negative.n.n.n.methylation<- n.n.n.methylation[, which(n.n.n.methylation[which(n.n.n.methylation$Probe.ID=="ER"), ]== "Negative")]
Negative.n.n.n.methylation$Probe.ID <- n.n.n.methylation$Probe.ID
Negative.n.n.n.methylation$Hugo_Symbol <- n.n.n.methylation$Hugo_Symbol
Negative.n.n.n.methylation$Entrez_Gene_Id <- n.n.n.methylation$Entrez_Gene_Id

Negative.n.n.n.methylation <- Negative.n.n.n.methylation %>% relocate(Entrez_Gene_Id, .before = TCGA.A2.A3XT.01)
Negative.n.n.n.methylation <- Negative.n.n.n.methylation %>% relocate(Hugo_Symbol, .before = Entrez_Gene_Id)
Negative.n.n.n.methylation <- Negative.n.n.n.methylation %>% relocate(Probe.ID, .before = Hugo_Symbol)

dim(Negative.n.n.n.methylation)

Negative.n.n.n.methylation[1:4, 1:5]

#======================================================================================================================================


save(Negative.n.n.n.data_CNA, file = "ER_test/Negative.n.n.n.data_CNA.2.RData")
save(Negative.n.n.n.expression, file = "ER_test/Negative.n.n.n.expression.2.RData")
save(Negative.n.n.n.methylation, file = "ER_test/Negative.n.n.n.methylation.2.RData")


n.n.data_CNA.2[1:4, 1:5]
n.n.expression.2[1:4, 1:5]
n.n.methylation.2[1:4, 1:5]

dim(n.n.data_CNA.2)
dim(n.n.expression.2)
dim(n.n.methylation.2)

#========================================================================================================================================

# Positive

#CNA
# Combine clinical data with main CNA
n.n.n.data_CNA <- rbind(n.n.data_CNA, t.n.ER.clinical.data)
dim(n.n.n.data_CNA)
n.n.n.data_CNA[1:5, 1:6]

# Identify ER column 
which(n.n.n.data_CNA$Probe.ID=="ER")

# Choose Positive 

Positive.n.n.n.data_CNA<- n.n.n.data_CNA[, which(n.n.n.data_CNA[which(n.n.n.data_CNA$Probe.ID=="ER"), ]== "Positive")]
Positive.n.n.n.data_CNA$Probe.ID <- n.n.n.data_CNA$Probe.ID
Positive.n.n.n.data_CNA$Hugo_Symbol <- n.n.n.data_CNA$Hugo_Symbol
Positive.n.n.n.data_CNA$Entrez_Gene_Id <- n.n.n.data_CNA$Entrez_Gene_Id

Positive.n.n.n.data_CNA <- Positive.n.n.n.data_CNA %>% relocate(Entrez_Gene_Id, .before = TCGA.A2.A0YK.01)
Positive.n.n.n.data_CNA <- Positive.n.n.n.data_CNA %>% relocate(Hugo_Symbol, .before = Entrez_Gene_Id)
Positive.n.n.n.data_CNA <- Positive.n.n.n.data_CNA %>% relocate(Probe.ID, .before = Hugo_Symbol)

dim(Positive.n.n.n.data_CNA)

Positive.n.n.n.data_CNA[1:4, 1:5]


#G.expression
# Combine clinical data with main CNA
n.n.n.expression <- rbind(n.n.expression, t.n.ER.clinical.data)
dim(n.n.n.expression)
n.n.n.expression[1:5, 1:6]

# Identify ER column 
which(n.n.n.expression$Probe.ID=="ER")

# Choose Positive 

Positive.n.n.n.expression<- n.n.n.expression[, which(n.n.n.expression[which(n.n.n.expression$Probe.ID=="ER"), ]== "Positive")]
Positive.n.n.n.expression$Probe.ID <- n.n.n.expression$Probe.ID
Positive.n.n.n.expression$Hugo_Symbol <- n.n.n.expression$Hugo_Symbol
Positive.n.n.n.expression$Entrez_Gene_Id <- n.n.n.expression$Entrez_Gene_Id

Positive.n.n.n.expression <- Positive.n.n.n.expression %>% relocate(Entrez_Gene_Id, .before = TCGA.A2.A0YK.01)
Positive.n.n.n.expression <- Positive.n.n.n.expression %>% relocate(Hugo_Symbol, .before = Entrez_Gene_Id)
Positive.n.n.n.expression <- Positive.n.n.n.expression %>% relocate(Probe.ID, .before = Hugo_Symbol)

dim(Positive.n.n.n.expression)

Positive.n.n.n.expression[1:4, 1:5]


#methylation
#Combine clinical data with main CNA
n.n.n.methylation<- rbind(n.n.methylation, t.n.ER.clinical.data)
dim(n.n.n.methylation)
n.n.n.methylation[1:5, 1:6]

# Identify ER column 
which(n.n.n.methylation$Probe.ID=="ER")

# Choose Positive 

Positive.n.n.n.methylation<- n.n.n.methylation[, which(n.n.n.methylation[which(n.n.n.methylation$Probe.ID=="ER"), ]== "Positive")]
Positive.n.n.n.methylation$Probe.ID <- n.n.n.methylation$Probe.ID
Positive.n.n.n.methylation$Hugo_Symbol <- n.n.n.methylation$Hugo_Symbol
Positive.n.n.n.methylation$Entrez_Gene_Id <- n.n.n.methylation$Entrez_Gene_Id

Positive.n.n.n.methylation <- Positive.n.n.n.methylation %>% relocate(Entrez_Gene_Id, .before = TCGA.A2.A0YK.01)
Positive.n.n.n.methylation <- Positive.n.n.n.methylation %>% relocate(Hugo_Symbol, .before = Entrez_Gene_Id)
Positive.n.n.n.methylation <- Positive.n.n.n.methylation %>% relocate(Probe.ID, .before = Hugo_Symbol)

dim(Positive.n.n.n.methylation)

Positive.n.n.n.methylation[1:4, 1:5]



save(Positive.n.n.n.data_CNA, file = "ER_test/Positive.n.n.n.data_CNA.2.RData")
save(Positive.n.n.n.expression, file = "ER_test/Positive.n.n.n.expression.2.RData")
save(Positive.n.n.n.methylation, file = "ER_test/Positive.n.n.n.methylation.2.RData")


n.n.data_CNA.2[1:4, 1:5]
n.n.expression.2[1:4, 1:5]
n.n.methylation.2[1:4, 1:5]

dim(n.n.data_CNA.2)
dim(n.n.expression.2)
dim(n.n.methylation.2)


