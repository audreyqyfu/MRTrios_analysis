#Title: Model identification by MRPC 

#Description: Using adjacency matrix to identify the different models of MRPC (i.e M0, M1...)  

#Date: June 1, 2021

#last updated: July 30, 2021


#=============================================================================================================

library(MRPC,lib="/mnt/ceph/megheib/Rpackages")


# Adj. matrix

full.adj.cor=read.delim("adj.cor.4.test.10.3.txt", sep = "\t", header = FALSE)
colnames(full.adj.cor) <- c("probe.id", "gene.name", "entrez.id",  "sub.ID" , "adj1", "adj2", "adj3", "adj4", "adj5", "adj6", "adj7", "adj8", "adj9","cor12", "cor13", "cor23")

dim(full.adj.cor)
full.adj.cor[1:10, ]


#================================================================
full.adj.cor.total.4 <- full.adj.cor

# apply MRPC with each trios

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

trios <- dim(full.adj.cor.total.4)[1]
print(trios)

n.data <- list()

corr <- list()

List.Match.significant.trios <- list()

id <- dim(full.adj.cor.total.4)[1]


# Construct trios 

for (i in 1:id) {
  
 #i=1
  
  full.adj.cor.total.4[i, ]
 
  S1 <-  full.adj.cor.total.4[i, ][, c(-1:-4, -14:-16)]
  
  full.adj.cor.total.4[1:5, c(-1:-4, -14:-16)]
  
  S2 <- matrix(unlist(S1), 3, 3)
  
  
  S2 <- as.data.frame(apply(S2, 2, as.numeric))
  
  S2 <- matrix(as.numeric(unlist(S2)), 3, 3)
  
 #print(S2)
  
  colnames(S2) <- c("V1", "T1","T2")
  
  
  Adj.infe1 <- S2
  
 # rownames(S2) <- rownames(Adj.M01)
  
 # print(S2)
  
  
  Adj.infe <- Adj.infe1[1:3,1:3] 
  
  rownames(Adj.infe) <- rownames(Adj.M01)
  
  #typeof(Adj.infe)
  
  #print(Adj.infe)
  
  colnames(Adj.infe) <- rownames(Adj.infe) <- colnames(Adj.M01) 
  
  if(identical(Adj.M0,Adj.infe) || identical(Adj.M01,Adj.infe)){
    M00.FDR.1[i]<-i
  }
  #print(Adj.M0)
  
 # print(M00.FDR.1[i])
  if(identical(Adj.M1,Adj.infe) || identical(Adj.M11,Adj.infe)){
    M11.FDR.1[i]<-i
  }
  
 # print(M11.FDR.1[i])
  ########
  
  if(identical(Adj.M1,Adj.infe) ){
    M11I.FDR.1I[i]<-i
  }
  
#  print( M11I.FDR.1I[i])
  
  if(identical(Adj.M11,Adj.infe)  ){
    M11.FDR.1II[i]<-i
  }
  
 # print(M11.FDR.1II[i])
  ###
  
  
  if(identical(Adj.M2,Adj.infe) || identical(Adj.M21,Adj.infe)){
    M22.FDR.1[i]<-i
  }
 # print( M22.FDR.1[i])
  
  if(identical(Adj.M3,Adj.infe)){
    M33.FDR.1[i]<-i
  }
 # print(M33.FDR.1[i])
  
  if(identical(Adj.M4,Adj.infe)){
    M44.FDR.1[i]<-i
  }
#  print(M44.FDR.1[i])
  
}



List.models.1.all.ADDIS <- list()
List.models.1.all.ADDIS$M0 <- M00.FDR.1[!is.na(M00.FDR.1)] 

List.models.1.all.ADDIS$M1 <- M11.FDR.1[!is.na(M11.FDR.1)]

List.models.1.all.ADDIS$type1 <- M11I.FDR.1I[!is.na(M11I.FDR.1I)]

List.models.1.all.ADDIS$type2 <- M11.FDR.1II[!is.na(M11.FDR.1II)]


List.models.1.all.ADDIS$M2 <- M22.FDR.1[!is.na(M22.FDR.1)] 
List.models.1.all.ADDIS$M3 <- M33.FDR.1[!is.na(M33.FDR.1)] 
List.models.1.all.ADDIS$M4 <- M44.FDR.1[!is.na(M44.FDR.1)] 


All.models.1.all.ADDIS.TCGA<- matrix(c(length(List.models.1.all.ADDIS$M0),
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
colnames(All.models.1.all.ADDIS.TCGA) <- c("M0","M1","M1I","M1II","M2","M3","M4","others", "probes")


print(All.models.1.all.ADDIS.TCGA)




write.table(All.models.1.all.ADDIS.TCGA, file = "All.models.1.all.ADDIS.TCGA.10.3.neg.txt")

read.delim("All.models.1.all.ADDIS.TCGA.10.txt")

