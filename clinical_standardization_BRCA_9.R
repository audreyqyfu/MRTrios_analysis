
#Title: Normalizing the methylation data 

# Description: Doing PCA normalization for methylation data, then applying MRPC to see the effect 

# Date: 01-01-2021

# Last updated: 03-12-2021

#==========================================================================================================
# Required package
require(ade4)
require(FactoMineR)
require(ggplot2)

# Reading CNA.exp.meth data 

n.data.CNA.exp.meth <- read.csv("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Cancer Genomic/brca_tcga/n.n.data.CNA.exp.meth.txt", sep="")
dim(n.data.CNA.exp.meth)
n.data.CNA.exp.meth[1:6, 1:6]
dim(n.data.CNA.exp.meth)

#extract meth data
data.Meth <- n.data.CNA.exp.meth[, seq(3, dim(n.data.CNA.exp.meth)[2],3) ]

data.Meth[1:10, 1:6]
dim(data.Meth)


#============================================================================================================

# Reading clinical data
clinical.data <- read.delim("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Cancer Genomic/brca_tcga/brca_tcga_clinical_data__grouping_5.txt", comment.char="#")
head(clinical.data)
dim(clinical.data)
names(clinical.data)

ER.clinical.data <- clinical.data[, -28]

# changing gene names, i.e TCGA-3C-AAAU-01 to TCGA.3C.AAAU.01

ER.clinical.data$Sample.ID <- gsub("\\-", ".", ER.clinical.data$Sample.ID)

head(ER.clinical.data)[1:3]

# matching clinical data with main data  

match.full.clinical <- match(rownames(n.data.CNA.exp.meth),ER.clinical.data$Sample.ID)

match.full.clinical.nona <-match.full.clinical[!is.na(match.full.clinical)]

length(match.full.clinical.nona)


n.ER.clinical.data <- ER.clinical.data[match.full.clinical.nona, ]

dim(n.ER.clinical.data)

dim(n.data.CNA.exp.meth)


#=======================================================================================================
# using ade4 to obtain PCs 

n.ER.clinical.data <- ER.clinical.data[match.full.clinical.nona, ]

dim(n.ER.clinical.data)


#naming the rows 

rownames(n.ER.clinical.data) <- n.ER.clinical.data[, 1] 

n.n.ER.clinical.data <- n.ER.clinical.data[, -1]

n.n.ER.clinical.data[1:5,1:3]

# Changing it to factors

#n.n.ER.clinical.data[] = apply(n.n.ER.clinical.data, 2, function(x) nlevels(as.factor(x)))

n.n.ER.clinical.data[] <- lapply(n.n.ER.clinical.data, factor)

# Using ade4

mca3 = dudi.acm(n.n.ER.clinical.data, scannf = FALSE, nf = 39)

mca3$eig

#PCA

dim(mca3$li)

head(mca3$li)

head(mca3$co)



# Save the PCA
PCA.clinical.data <- mca3$li
save(PCA.clinical.data,file="PCA.clinical.data.RData")
write.csv(PCA.clinical.data,file="PCA.clinical.data.csv")

#==========================================================================================================
# Normalizing Meth data 

# Reading 

#load("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Cancer Genomic/MCA_PCA/PCA.clinical.data.RData")


head(PCA.clinical.data)[1:5]

dim(PCA.clinical.data)

PCA.stand.meth <- matrix(0,dim(data.Meth)[1], dim(data.Meth)[2] )

PCA.stand.meth.1 <- matrix(0,dim(data.Meth)[1], 9 )

dim(PCA.stand.meth)

## Residuals of PCA


for (i in 1:9) {
  
  model <- lm(log(data.Meth[,i]/(1-data.Meth[,i]))~., data =PCA.clinical.data )
  summary(model)
  PCA.stand.meth.1[, i] =model$residuals
  
  # print( PCA.stand.meth[,i])
}



# naming 

colnames(PCA.stand.meth) <- colnames(data.Meth)
rownames(PCA.stand.meth) <- rownames(data.Meth)


data.Meth[1:3, 1:6]

PCA.stand.meth.1[1:3, 1:6]


#Plotting normalized ones 

pairs(PCA.stand.meth[,1:9], pch=16, cex=0.6, col="blue")

#Original data 

pairs(data.Meth[,1:9],pch=16, cex=0.6, col="black" )


#===========================================================================================================
#Replacing the original variable by standardized ones


n.data.CNA.exp.peer.meth <- replace(n.data.CNA.exp.meth, seq(2, dim(n.data.CNA.exp.meth)[2],3), peer.expression )

n.data.CNA.exp.peer.meth[1:6, 1:6]

n.data.CNA.exp.peer.meth.PCA <- replace(n.data.CNA.exp.peer.meth, seq(3, dim(n.data.CNA.exp.peer.meth)[2],3), PCA.stand.meth )

n.data.CNA.exp.peer.meth.PCA[1:6, 1:6]

dim(n.data.CNA.exp.peer.meth.PCA)

############################################Applying MRPC on the new data sets########################################################

#General results


library(MRPC)

# apply MRPC with adjusted PCs on each trios

#truth for M
#V1-->T1
Truth.M0 <- MRPCtruth$M0
Adj.M0<- as(Truth.M0,"matrix")
#V1-->T2
Adj.M01 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M01) <- colnames(Adj.M01) <- colnames(Adj.M0)
Adj.M01[1,3] <- 1

#truth for M1
#V1-->T1-->T2
Truth.M1 <- MRPCtruth$M1
Adj.M1<- as(Truth.M1,"matrix")
#V1-->T2-->T1
Adj.M11 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M11) <- colnames(Adj.M11) <- colnames(Adj.M1)
Adj.M11[1,3] <- 1
Adj.M11[3,2] <- 1

#truth for M2
#V1-->T1<--T2
Truth.M2 <- MRPCtruth$M2
Adj.M2<- as(Truth.M2,"matrix")
#V1-->T2<--T1
Adj.M21 <- matrix(0,nrow=3,ncol = 3)
rownames(Adj.M21) <-colnames(Adj.M21) <-colnames(Adj.M2)
Adj.M21[1,3] <- 1
Adj.M21[2,3] <- 1
#truth for M3
#V1-->T1, V1-->T2
Truth.M3 <- MRPCtruth$M3
Adj.M3 <- as(Truth.M3,"matrix")

#truth for M4
#V1-->T1, V1-->T2, T1--T2
Truth.M4 <- MRPCtruth$M4
Adj.M4 <- as(Truth.M4,"matrix")


#Initial results
M00.FDR.1 <- c()
M11.FDR.1 <- c()
M11I.FDR.1I <- c()
M11.FDR.1II <- c()

M22.FDR.1 <- c()
M33.FDR.1 <- c()
M44.FDR.1 <- c()
MRPC.fit.table<- list()
List.Associated.PCs.1 <- list()

# columns number for trios
start.col <- seq(1,dim(n.data.CNA.exp.peer.meth.PCA)[2],3)
end.col <- seq(3,dim(n.data.CNA.exp.peer.meth.PCA)[2],3)

trios <- dim(n.data.CNA.exp.peer.meth.PCA)[2]/3
print(trios)



id <- length(n.data.CNA.exp.peer.meth.PCA[,1])

n.data <- list()

List.Match.significant.trios <- list()

for (ii in 1:id) {
  # data
  
  data <- n.data.CNA.exp.peer.meth.PCA[,start.col[ii]:end.col[ii]]
  # search in each PCs
  head(data)
  #data <- simu_data_M0
  
  n <- nrow (data)
  
  
  colnames(data) <- c("V1", "T1","T2")
  
  V <- colnames(data)     # Column names
  
  
  # Classical correlation
  suffStat <- list(C = cor(data,use = "complete.obs"),
                   n = n)
  
  MRPC.fit.FDR.ADDIS <- MRPC(data,
                             suffStat,
                             GV = 1,
                             FDR = 0.05,
                             indepTest = 'gaussCItest',
                             labels = V,
                             FDRcontrol = "ADDIS",
                             verbose = FALSE)
  
  
  #plot(MRPC.fit_FDR)
  Adj.infe1 <- as( MRPC.fit.FDR.ADDIS@graph,"matrix")
  Adj.infe <- Adj.infe1[1:3,1:3] #only consider snp, cis, trans
  colnames(Adj.infe) <- rownames(Adj.infe) <- colnames(Adj.M01) 
  
  if(identical(Adj.M0,Adj.infe) || identical(Adj.M01,Adj.infe)){
    M00.FDR.1[ii]<-ii
  }
  
  if(identical(Adj.M1,Adj.infe) || identical(Adj.M11,Adj.infe)){
    M11.FDR.1[ii]<-ii
  }
  
  ########
  
  if(identical(Adj.M1,Adj.infe) ){
    M11I.FDR.1I[ii]<-ii
  }
  
  if(identical(Adj.M11,Adj.infe)  ){
    M11.FDR.1II[ii]<-ii
  }
  
  
  ###
  
  
  if(identical(Adj.M2,Adj.infe) || identical(Adj.M21,Adj.infe)){
    M22.FDR.1[ii]<-ii
  }
  if(identical(Adj.M3,Adj.infe)){
    M33.FDR.1[ii]<-ii
  }
  if(identical(Adj.M4,Adj.infe)){
    M44.FDR.1[ii]<-ii
  }
  
}



List.models.1.all.ADDIS <- list()
List.models.1.all.ADDIS$M0 <- M00.FDR.1[!is.na(M00.FDR.1)] 

List.models.1.all.ADDIS$M1 <- M11.FDR.1[!is.na(M11.FDR.1)]

List.models.1.all.ADDIS$type1 <- M11I.FDR.1I[!is.na(M11I.FDR.1I)]

List.models.1.all.ADDIS$type2 <- M11.FDR.1II[!is.na(M11.FDR.1II)]


List.models.1.all.ADDIS$M2 <- M22.FDR.1[!is.na(M22.FDR.1)] 
List.models.1.all.ADDIS$M3 <- M33.FDR.1[!is.na(M33.FDR.1)] 
List.models.1.all.ADDIS$M4 <- M44.FDR.1[!is.na(M44.FDR.1)] 
#List.models.1.all.ADDIS$M1 <- MRPC.fit.table[!is.na(MRPC.fit.table)] 


All.models.1.all.ADDIS<- matrix(c(length(List.models.1.all.ADDIS$M0),
                                  length(List.models.1.all.ADDIS$M1),
                                  length(List.models.1.all.ADDIS$type1),
                                  length(List.models.1.all.ADDIS$type2),
                                  length(List.models.1.all.ADDIS$M2),
                                  length(List.models.1.all.ADDIS$M3),
                                  length(List.models.1.all.ADDIS$M4), 
                                  length(setdiff(1:trios,c(List.models.1.all.ADDIS$M0,List.models.1.all.ADDIS$M1,
                                                           List.models.1.all.ADDIS$M2,List.models.1.all.ADDIS$M3,
                                                           List.models.1.all.ADDIS$M4))), 
                                  
                                  id  ),
                                nrow = 1,ncol = 9)
colnames(All.models.1.all.ADDIS) <- c("M0","M1","M1I","M1II","M2","M3","M4","others", "Gene")


print(All.models.1.all.ADDIS)







