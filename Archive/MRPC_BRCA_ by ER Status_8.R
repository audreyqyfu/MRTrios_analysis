
#Title: Applying MRPC on the BRCA data by the ER Status 

# Description: Dividing the data by ER test status then applying MRPC 

#Date: 02-20-2021

#Last updated: 04-25-2021


library(MRPC)

# Reading the data

data_CNA <- read.delim("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Cancer Genomic/brca_tcga/data_CNA.txt")
head(data_CNA)[1:10]
dim(data_CNA)


methylation <- read.delim("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Cancer Genomic/brca_tcga/data_methylation_hm450.txt")
head(methylation)[1:10]
dim(methylation)


expression <- read.delim("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Cancer Genomic/brca_tcga/data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt")
head(expression)[1:10]
dim(expression)



#CNA/methylation_rows
dim(data_CNA)
dim(methylation)
dim(expression)

#meth/CNA
match.meth.CNA <- match(methylation$Hugo_Symbol, data_CNA$Hugo_Symbol)

length(match.meth.CNA)

match.meth.CNA.nona <- match.meth.CNA[!is.na(match.meth.CNA)]

length(match.meth.CNA.nona)

n.data_CNA <- data_CNA[match.meth.CNA.nona, ]

dim(n.data_CNA)

#CNA/meth
match.CNA.meth <- match(data_CNA$Hugo_Symbol, methylation$Hugo_Symbol)

length(match.CNA.meth)

match.CNA.meth.nona <- match.CNA.meth [!is.na(match.CNA.meth )]

length(match.CNA.meth.nona)

n.methylation <- methylation[match.CNA.meth.nona, ]

dim(n.methylation)


## CNA/methylation_Match column

#meth/CNA
match.indv.meth.CNA <- match(colnames(methylation),colnames(data_CNA))

match.indv.meth.CNA.nona <- match.indv.meth.CNA[!is.na(match.indv.meth.CNA)]

length(match.indv.meth.CNA.nona)

#CNA/meth
match.indv.CNA.meth<- match(colnames(data_CNA),colnames(methylation))

match.indv.CNA.meth.nona <- match.indv.CNA.meth[!is.na(match.indv.CNA.meth)]

length(match.indv.CNA.meth.nona)


# new matched data sets

n.n.data_CNA <- data_CNA[match.meth.CNA.nona, match.indv.meth.CNA.nona]

dim(n.n.data_CNA)

head(n.n.data_CNA)[1:6]

#head(n.n.data_CNA[order(n.n.data_CNA$Hugo_Symbol),])[1:6]


n.n.methylation <- methylation[match.CNA.meth.nona,match.indv.CNA.meth.nona ]

dim(n.n.methylation)

head(n.n.methylation)[1:6]

#head(n.n.methylation[order(n.n.methylation$Hugo_Symbo),])[1:6]

#Checking the dimensions

names(n.n.methylation)==names(n.n.data_CNA)

sort(n.n.methylation$Hugo_Symbol)==sort(n.n.data_CNA$Hugo_Symbol)


#############################################3#expression/methylation_rows###############################

dim(n.n.data_CNA)
dim(n.n.methylation)

dim(expression)
dim(methylation)

#meth/expression
match.meth.exp <- match(n.n.methylation$Hugo_Symbol, expression$Hugo_Symbol)

length(match.meth.exp)

match.meth.exp.nona <- match.meth.exp[!is.na(match.meth.exp)]

length(match.meth.exp.nona)


n.expression <- expression[match.meth.exp.nona, ]
  
dim(n.expression)


#expression/meth
match.exp.meth <- match(n.expression$Hugo_Symbol, n.n.methylation$Hugo_Symbol)

length(match.exp.meth)

match.exp.meth.nona <- match.exp.meth[!is.na(match.exp.meth)]

length(match.exp.meth.nona)


# expression/methylation_Match column

match.indv.meth.exp <- match(colnames(n.n.methylation),colnames(expression))

match.indv.meth.exp.nona <- match.indv.meth.exp[!is.na(match.indv.meth.exp)]

length(match.indv.meth.exp.nona)


# new matched data sets of expression 
n.n.expression.1 <- expression[match.exp.meth.nona,match.indv.CNA.meth.nona]
dim(n.n.expression.1)
n.n.expression<- n.n.expression.1[!duplicated(n.n.expression.1[, c("Hugo_Symbol")]),]
dim(n.n.expression)


#=========================================================================================================
#dim(n.data_CNA)
#dim(n.methylation)
#dim(n.expression)


dim(n.n.data_CNA)
dim(n.n.methylation)
dim(n.n.expression)

#n.n.methylation[order(n.n.methylation$Hugo_Symbol),][1:10,1:3]
#n.n.data_CNA[order(n.n.data_CNA$Hugo_Symbol),][1:10,1:3]
#n.n.expression[order(n.n.expression$Hugo_Symbol),][1:10,1:3]


# intersection 

id <- intersect(unique(n.n.expression$Hugo_Symbol),n.n.methylation$Hugo_Symbol)
length(id)

sub.expression <- subset(n.n.expression,unique(n.n.expression$Hugo_Symbol) %in% unique(id))

dim(sub.expression)

sub.data_CNA <- subset(n.n.data_CNA,n.n.data_CNA$Hugo_Symbol %in% id)

dim(sub.data_CNA)


sub.methylation <- subset(n.n.methylation,n.n.methylation$Hugo_Symbol %in% id)

dim(sub.methylation)



sub.methylation[order(sub.methylation$Hugo_Symbol),][1:10,1:3]
sub.data_CNA[order(sub.data_CNA$Hugo_Symbol),][1:10,1:3]
sub.expression[order(sub.expression$Hugo_Symbol),][1:10,1:3]


n.n.methylation <- sub.methylation
n.n.data_CNA <- sub.data_CNA
n.n.expression <- sub.expression

dim(n.n.data_CNA)
dim(n.n.methylation)
dim(n.n.expression)

head(n.n.data_CNA[, 1])
length(n.n.data_CNA[, 1])

# comparing
comparison <- compare(n.n.data_CNA[, 1],n.n.methylation[, 1], allowAll=TRUE)
length(comparison$tM)
comparison <- compare(n.n.data_CNA[, 1],n.n.expression[, 1], allowAll=TRUE)
length(comparison$tM)
comparison <- compare(n.n.methylation[, 1],n.n.expression[, 1], allowAll=TRUE)
length(comparison$tM)

#===========================================================================================================
# Constructing trios 

full.data <- c(n.n.data_CNA[,1],n.n.expression[,1],n.n.methylation[,1])

full.data <- n.n.methylation[,1]
head(full.data)

dim(full.data)


id <- length(full.data)

CNA.exp.meth <- NULL
CNA <- NULL
exp <- NULL
meth <- NULL

for (i in 1:id) {
 # i=2
  CNA[i] <- paste0(full.data[i])
  exp[i] <- paste0(full.data[i])
  meth[i] <- paste0(full.data[i])
  
  CNA.exp.meth[[i]]<-c(n.n.data_CNA[which(n.n.data_CNA$Hugo_Symbol==CNA[i])[1],],
                       n.n.expression[which(n.n.expression$Hugo_Symbol==exp[i])[1],],
                       n.n.methylation[which(n.n.methylation$Hugo_Symbol==meth[i])[1],])
}

Col.names <- c(CNA, exp, meth)
data.CNA.exp.meth <- matrix(as.numeric(unlist(CNA.exp.meth)),nrow=dim(n.n.expression)[2],ncol=id*3, byrow = FALSE)

seq.list <- NULL
for(j in 1:id){
  seq.list[[j]]<-c(seq(j, id*3, id))
}
colnames(data.CNA.exp.meth) <- Col.names[as.numeric(unlist(seq.list))]
rownames(data.CNA.exp.meth) <- colnames(n.n.expression)
dim(data.CNA.exp.meth)
#data.CNA.exp.meth[3:10,1:10]

n.data.CNA.exp.meth <- data.CNA.exp.meth[-1:-2,]

n.data.CNA.exp.meth[1:6, 1:6]

dim(n.data.CNA.exp.meth)



######################Adding demographics to the main data###############################################


# Reading the clinical data 

clinical.data <- read.delim("C:/Users/molot/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Cancer Genomic/brca_tcga/brca_tcga_clinical_data.txt", comment.char="#")

head(clinical.data)
dim(clinical.data)
names(clinical.data)


ER.clinical.data <- clinical.data[, c(3,28)]

dim(ER.clinical.data)

# changing gene names, i.e TCGA-3C-AAAU-01 to TCGA.3C.AAAU.01

ER.clinical.data$Sample.ID <- gsub("\\-", ".", ER.clinical.data$Sample.ID)

head(ER.clinical.data)


match.full.clinical <- match(rownames(n.data.CNA.exp.meth),ER.clinical.data$Sample.ID)



match.full.clinical.nona <-match.full.clinical[!is.na(match.full.clinical)]

length(match.full.clinical.nona)


n.ER.clinical.data <- ER.clinical.data[match.full.clinical.nona, ]

dim(n.ER.clinical.data)

length(n.ER.clinical.data$Sample.ID)==length(rownames(n.data.CNA.exp.meth))

# Changing data frame

n.data.CNA.exp.meth.frame <- data.frame(n.data.CNA.exp.meth)

n.data.CNA.exp.meth.frame$Sample.ID <- n.ER.clinical.data$Sample.ID


# merging the data
n.n.data.CNA.exp.meth.full <- merge(ER.clinical.data, n.data.CNA.exp.meth.frame, by = "Sample.ID" )


#n.n.data.CNA.exp.meth[which(n.n.data.CNA.exp.meth$ER.Status.By.IHC== "Negative"),][1:6,1:8]


Negative.n.n.data.CNA.exp.meth <- n.n.data.CNA.exp.meth.full[which(n.n.data.CNA.exp.meth.full$ER.Status.By.IHC== "Negative"),]

Positive.n.n.data.CNA.exp.meth <- n.n.data.CNA.exp.meth.full[which(n.n.data.CNA.exp.meth.full$ER.Status.By.IHC== "Positive"),]

n.n.data.CNA.exp.meth.full[1:6,1:8]

Negative.n.n.data.CNA.exp.meth[1:6,1:8]

Positive.n.n.data.CNA.exp.meth[1:6,1:8]

n.n.data.CNA.exp.meth <- n.n.data.CNA.exp.meth.full[, -1:-2]

###############################General results ###########################################################################


###############################General results #######

n.n.data.CNA.exp.meth[1:6,1:8]



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
start.col <- seq(1,dim(n.n.data.CNA.exp.meth)[2],3)
end.col <- seq(3,dim(n.n.data.CNA.exp.meth)[2],3)

trios <- dim(n.n.data.CNA.exp.meth)[2]/3
print(trios)


n.data <- list()

List.Match.significant.trios <- list()

for (ii in 1:id) {
  # data


  data <- n.n.data.CNA.exp.meth[,start.col[ii]:end.col[ii]]
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



############################Negative#########################################################################################


#Negative.n.n.data.CNA.exp.meth <- n.n.data.CNA.exp.meth[which(n.n.data.CNA.exp.meth$ER.Status.By.IHC== "Negative"),]


Negative.n.n.data.CNA.exp.meth[1:10,1:10]


n.Negative.n.n.data.CNA.exp.meth <- Negative.n.n.data.CNA.exp.meth[,-1:-2]

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
start.col <- seq(1,dim(n.Negative.n.n.data.CNA.exp.meth)[2],3)
end.col <- seq(3,dim(n.Negative.n.n.data.CNA.exp.meth)[2],3)

trios <- dim(n.Negative.n.n.data.CNA.exp.meth)[2]/3
print(trios)


n.data <- list()

List.Match.significant.trios <- list()

for (ii in 1:id) {
  # data
  
  data <- n.Negative.n.n.data.CNA.exp.meth[,start.col[ii]:end.col[ii]]
  # search in each PCs
  
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
  
  # plot(MRPC.fit.FDR.ADDIS)
  #print(MRPC.fit.FDR.ADDIS)
  
  
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



List.models.1.all.ADDIS.Negative <- list()
List.models.1.all.ADDIS.Negative$M0 <- M00.FDR.1[!is.na(M00.FDR.1)] 

List.models.1.all.ADDIS.Negative$M1 <- M11.FDR.1[!is.na(M11.FDR.1)]

List.models.1.all.ADDIS.Negative$type1 <- M11I.FDR.1I[!is.na(M11I.FDR.1I)]

List.models.1.all.ADDIS.Negative$type2 <- M11.FDR.1II[!is.na(M11.FDR.1II)]


List.models.1.all.ADDIS.Negative$M2 <- M22.FDR.1[!is.na(M22.FDR.1)] 
List.models.1.all.ADDIS.Negative$M3 <- M33.FDR.1[!is.na(M33.FDR.1)] 
List.models.1.all.ADDIS.Negative$M4 <- M44.FDR.1[!is.na(M44.FDR.1)] 
#List.models.1.all.ADDIS.Negative$M1 <- MRPC.fit.table[!is.na(MRPC.fit.table)] 


All.models.1.all.ADDIS.Negative<- matrix(c(length(List.models.1.all.ADDIS.Negative$M0),
                                     length(List.models.1.all.ADDIS.Negative$M1),
                                     length(List.models.1.all.ADDIS.Negative$type1),
                                     length(List.models.1.all.ADDIS.Negative$type2),
                                     length(List.models.1.all.ADDIS.Negative$M2),
                                     length(List.models.1.all.ADDIS.Negative$M3),
                                     length(List.models.1.all.ADDIS.Negative$M4), 
                                     length(setdiff(1:trios,c(List.models.1.all.ADDIS.Negative$M0,List.models.1.all.ADDIS.Negative$M1,
                                                              List.models.1.all.ADDIS.Negative$M2,List.models.1.all.ADDIS.Negative$M3,
                                                              List.models.1.all.ADDIS.Negative$M4))), 
                                     
                                     id  ),
                                   nrow = 1,ncol = 9)

colnames(All.models.1.all.ADDIS.Negative) <- c("M0","M1","M1I","M1II","M2","M3","M4","others", "Gene")

print(All.models.1.all.ADDIS.Negative)



#############################################Positive#############################################################################


#Positive.n.n.data.CNA.exp.meth <- n.n.data.CNA.exp.meth[which(n.n.data.CNA.exp.meth$ER.Status.By.IHC== "Positive"),]



Positive.n.n.data.CNA.exp.meth[1:10,1:10]

n.Positive.n.n.data.CNA.exp.meth <- Positive.n.n.data.CNA.exp.meth[,-1:-2]




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
start.col <- seq(1,dim(n.Positive.n.n.data.CNA.exp.meth)[2],3)
end.col <- seq(3,dim(n.Positive.n.n.data.CNA.exp.meth)[2],3)

trios <- dim(n.Positive.n.n.data.CNA.exp.meth)[2]/3
print(trios)


n.data <- list()

List.Match.significant.trios <- list()

for (ii in 1:id) {
  # data

  data <- n.Positive.n.n.data.CNA.exp.meth[,start.col[ii]:end.col[ii]]
  # search in each PCs
  
  #data <- simu_data_M0
  
  n <- nrow (data)
  
  
  colnames(data) <- c("V1", "T1","T2")
  
  V <- colnames(data)     # Column names
  
  
  
  #n.data[[ii]] <- data
  
  # data <- n.data[[ii]]
  
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
  
  MRPC.fit.table[[ii]]=as(MRPC.fit.FDR.ADDIS@graph, "matrix")
  
  
  
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



List.models.1.all.ADDIS.Positive <- list()
List.models.1.all.ADDIS.Positive$M0 <- M00.FDR.1[!is.na(M00.FDR.1)] 

List.models.1.all.ADDIS.Positive$M1 <- M11.FDR.1[!is.na(M11.FDR.1)]

List.models.1.all.ADDIS.Positive$type1 <- M11I.FDR.1I[!is.na(M11I.FDR.1I)]

List.models.1.all.ADDIS.Positive$type2 <- M11.FDR.1II[!is.na(M11.FDR.1II)]


List.models.1.all.ADDIS.Positive$M2 <- M22.FDR.1[!is.na(M22.FDR.1)] 
List.models.1.all.ADDIS.Positive$M3 <- M33.FDR.1[!is.na(M33.FDR.1)] 
List.models.1.all.ADDIS.Positive$M4 <- M44.FDR.1[!is.na(M44.FDR.1)] 
#List.models.1.all.ADDIS.Positive$M1 <- MRPC.fit.table[!is.na(MRPC.fit.table)] 


All.models.1.all.ADDIS.Positive<- matrix(c(length(List.models.1.all.ADDIS.Positive$M0),
                                     length(List.models.1.all.ADDIS.Positive$M1),
                                     length(List.models.1.all.ADDIS.Positive$type1),
                                     length(List.models.1.all.ADDIS.Positive$type2),
                                     length(List.models.1.all.ADDIS.Positive$M2),
                                     length(List.models.1.all.ADDIS.Positive$M3),
                                     length(List.models.1.all.ADDIS.Positive$M4), 
                                     length(setdiff(1:trios,c(List.models.1.all.ADDIS.Positive$M0,List.models.1.all.ADDIS.Positive$M1,
                                                              List.models.1.all.ADDIS.Positive$M2,List.models.1.all.ADDIS.Positive$M3,
                                                              List.models.1.all.ADDIS.Positive$M4))), 
                                     
                                     id  ),
                                   nrow = 1,ncol = 9)
colnames(All.models.1.all.ADDIS.Positive) <- c("M0","M1","M1I","M1II","M2","M3","M4","others", "Gene")

All.models.1.all.ADDIS.Positive
#####################





#Overall results
print(All.models.1.all.ADDIS)

# Negative results
print(All.models.1.all.ADDIS.Negative)


# Positive results
print(All.models.1.all.ADDIS.Positive)





