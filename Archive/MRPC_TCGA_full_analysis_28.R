#Title: TCGA analysis by MRPC 

#Description: We do cleaning and preparation for the data, then apply MRPC 

#created by Mohamed Megheib

#Date: 07-01-2021

#last updated: 09-01-2021


#=========================================================================================================================

library(MRPC,lib="/mnt/ceph/megheib/Rpackages")

#require(dplyr)

#=========================================================================================================================


# Reading CNA,  Gene expression and methylation  

load("CNA.t.10.unique.RData")
load("EM.t.10.RData")
load("MM.t.10.RData")

#Checking dimensions

dim(CNA.t.10.unique)
dim(EM.t.10)
dim(MM.t.10)

# Deleting rows with all-NAs from meythlation 

MM.t.10= MM.t.10[-which(rowSums(is.na(MM.t.10[ , 4:773]))==770), ]

dim(MM.t.10)
MM.t.10[1:5, 1:6]


# Deleting rows with all-NAs from gene.expression 

EM.t.10= EM.t.10[-which(rowSums(is.na(EM.t.10[ , 4:773]))==770), ]
dim(EM.t.10)
EM.t.10[1:5, 1:6]


#====================================================================================================================================
#choose non-idential gene expression indiviudals - Exclusing those gene whose maximum frequency >500

expression <- EM.t.10

freq <- matrix(0, dim(expression)[1], 1)

for (i in 1:dim(expression)[1]) {
  
  if(max(table(as.numeric(expression[i,-1:-3])))>500) {freq[i]=1}
  
  else{freq[i]=0}
  
}

sum(freq)


max(table(as.numeric(expression[i,-1:-3])))>500

#Deleting idential individuals

expression.2 <- expression[-which(freq %in% c(1)), ]
dim(expression.2)

EM.t.10.nonidentical <- expression.2

#save(EM.t.10.nonidentical, file = "EM.t.10.nonidentical.RData")

#================================================================================================================================
# Applying logit function on meythlation

my_fun <- function(x){log(x/(1-x))}

# Apply function to specific columns
data_apply <- apply(MM.t.10[ , 4:773 ], 2, my_fun)   

data_apply[1:4, 1:5] 

# Replicate original data
data_new <- MM.t.10        

# Replace specific columns
data_new[ , colnames(data_new) %in% colnames(data_apply)] <- data_apply  
data_new[1:4, 1:5] 

#===================================================================================================================

n.n.data_CNA.2 <- CNA.t.10.unique
n.n.expression.2 <- expression.2
n.n.methylation.3 <- data_new
dim(n.n.methylation.3)

# Search for NAs entrez ids in CNA data

full.data.names <- n.n.methylation.3$Hugo_Symbol
id <- length(full.data.names)

CNA.list <- n.n.data_CNA.2[,2:3]

for (i in 1:id) {
  
  if (is.na(n.n.methylation.3$Entrez_Gene_Id[i])==TRUE){
    
    meth <- NULL
    
    meth[i] <- paste0(full.data.names[i])
    n.n.methylation.3$Entrez_Gene_Id[i] <- matrix(c(unlist((CNA.list[which(CNA.list$Hugo_Symbol==meth[i]),]))), 1, 2)[, 2]
    
    }
}

# Search for NAs entrez ids in gene expression data

full.data.names <- n.n.methylation.3$Hugo_Symbol
expression.list <- n.n.expression.2[,2:3]
head(expression.list)

for (i in 1:id) {
  
  if (is.na(n.n.methylation.3$Entrez_Gene_Id[i])==TRUE){
    
    meth <- NULL
    
    meth[i] <- paste0(full.data.names[i])
    
    n.n.methylation.3$Entrez_Gene_Id[i] <- matrix(c(unlist((expression.list[which(expression.list$Hugo_Symbol==meth[i]),]))), 1, 2)[, 2]

  }
}



n.n.methylation.3.na <- n.n.methylation.3[is.na(n.n.methylation.3$Entrez_Gene_Id), ]
dim(n.n.methylation.3.na)

# Search for NAs entrez ids in master entrez id list 

my.symbols <- n.n.methylation.3$Hugo_Symbol

load("ntrez_id.2.RData")

full.data.names <- n.n.methylation.3$Hugo_Symbol
ntrez_id.list <- unique(ntrez_id.2)
head(ntrez_id.list)

for (i in 1:id) {
  
  if (is.na(n.n.methylation.3$Entrez_Gene_Id[i])==TRUE){
    
    meth <- NULL
    
    meth[i] <- paste0(full.data.names[i])
    
    n.n.methylation.3$Entrez_Gene_Id[i] <- matrix(c(unlist((ntrez_id.list[which(ntrez_id.list$Hugo_Symbol==meth[i]),]))), 1, 2)[, 2]
    
  }
}


n.n.methylation.3.na <- n.n.methylation.3[is.na(n.n.methylation.3$Entrez_Gene_Id), ]
dim(n.n.methylation.3.na)

n.n.data_CNA.2[1:4, 1:5]
#sum(is.na(n.n.data_CNA.2$Entrez_Gene_Id))
n.n.expression.2[1:4, 1:5]
#sum(is.na(n.n.expression.2$Entrez_Gene_Id))
n.n.methylation.3[1:4, 1:5]

n.n.methylation.3.na[1:4, 1:5]



dim(n.n.data_CNA.2)
dim(n.n.expression.2)
dim(n.n.methylation.3)
dim(n.n.methylation.3.na)

save(n.n.data_CNA.2, file = "TCGA_PCA/n.n.data_CNA.4.final.RData")
save(n.n.expression.2, file = "TCGA_PCA/n.n.expression.4.final.RData")
save(n.n.methylation.3, file = "TCGA_PCA/n.n.methylation.4.final.RData")

#=============================================================================================================================================    

# Creating a function to match trios and apply MRPC 

# Defining new variables 

full.range <- unique(n.n.methylation.3$Hugo_Symbol)
full.range.2 <- unique(n.n.methylation.3$Entrez_Gene_Id)


length(full.range.2)
full.range.2 [20:22]

which(is.na(full.range.2))

f<-function(full.range) {
  
  # when entriz id is NA, we match based on gene names
  
if (is.na(full.range)==TRUE){
    
full.range <- unique(n.n.methylation.3.na$Hugo_Symbol)
    
id <- length(full.range)
    
telliter <- 2000
    
#Construct trios and apply MRPC
    
for (i in 1:id) {
      
# showing number of iterations every telliter
    
if( i %% telliter == 0 ) cat(paste("iteration", i, "complete\n"))
      
CNA.exp.meth <- matrix(0)
CNA <- NULL
exp <- NULL
meth <- NULL
      
CNA[i] <- paste0(full.range[i])
      
exp[i] <- paste0(full.range[i])
meth[i] <- paste0(full.range[i])
      
# Adding filters to skip when the data do not match
      
if(dim(n.n.data_CNA.2[which(n.n.data_CNA.2$Hugo_Symbol==CNA[i]),])[1]==0) next
if(dim(n.n.expression.2[which(n.n.expression.2$Hugo_Symbol==exp[i]),])[1]==0) next
      
# defining matrices to save results 
      
total.probes <- list()

adj.total <- matrix(0, dim(n.n.methylation.3.na[which(n.n.methylation.3.na$Hugo_Symbol==meth[i]),])[1], 9)
      
adj.cor <- matrix(0, dim(n.n.methylation.3.na[which(n.n.methylation.3.na$Hugo_Symbol==meth[i]),])[1], 16)
      
# Creating a loop to include repeated cases 
      
for (j in 1: dim(n.n.methylation.3.na[which(n.n.methylation.3.na$Hugo_Symbol==meth[i]),])[1]) {

if(dim(n.n.methylation.3.na[which(n.n.methylation.3.na$Hugo_Symbol==meth[i]),][j, ])[1]==0) next
        
        
CNA.exp.meth<-matrix(c(n.n.data_CNA.2[which(n.n.data_CNA.2$Hugo_Symbol==CNA[i]),],
                               n.n.expression.2[which(n.n.expression.2$Hugo_Symbol==exp[i]),],
                               n.n.methylation.3.na[which(n.n.methylation.3.na$Hugo_Symbol==meth[i]),][j, ]), dim(n.n.expression.2)[2] , 3)
        
 #print(CNA.exp.meth[1:6, 1:3])
   
total.probes[j] <- CNA.exp.meth[1,3]
        
        
CNA.exp.meth[1:4, 1:3]
        
CNA.exp.meth.1 <- CNA.exp.meth[-1:-3,]
        
CNA.exp.meth.1[1:4, 1:3]  
        
if(sum(is.na(CNA.exp.meth.1[, 3]))==dim(CNA.exp.meth.1)[1]) next 
        
CNA.exp.meth.2 <- matrix(as.numeric(unlist(CNA.exp.meth.1)), dim(CNA.exp.meth.1)[1], 3)
        
CNA.exp.meth.2[1:4, 1:3]
data <-CNA.exp.meth.2
        

 dim(data)
        
 n <- nrow (data)
        
        
colnames(data) <- c("V1", "T1","T2")
        
 V <- colnames(data)     # Column names
        
#  print(cor(data, use = "complete.obs"))
        
cor12 <- cor(data[, 1], data[, 2], use = "complete.obs")
cor13 <- cor(data[, 1], data[, 3], use = "complete.obs")
cor23 <- cor(data[, 2], data[, 3], use = "complete.obs")
        
        
        
# Classical correlation
suffStat <- list(C = cor(data,use = "complete.obs"),
                         n = n)
        
# Corrected correlation
suffStat$C[1,2] <- suffStat$C[2,1] <- cor12 
suffStat$C[1,3] <- suffStat$C[3,1] <- cor13 
suffStat$C[2,3] <- suffStat$C[3,2] <- cor23
        
        
        
cor.mtx <- suffStat$C
        
        
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
adj.total[j, ] <- matrix(Adj.infe1, 1, 9)
        
        
adj.cor[j, ] <- c(unlist(CNA.exp.meth[1, 3]), unlist(CNA.exp.meth[2, 3]),unlist(CNA.exp.meth[3, 3]), j,adj.total[j,], cor.mtx[1,2], cor.mtx[1,3], cor.mtx[2,3])
        

#  print(t(adj.cor[j, ]))

 write.table(t(adj.cor[j, ]), file="adj.cor.4.test.10.31.txt", sep = "\t", row.names = FALSE,col.names = FALSE , append =TRUE,quote=FALSE)
        
             
      }
    }
    
  } else{
    

 #when entriz id is not NA, we match based on it 

    
id <- length(full.range)
    
telliter <- 2001

#Construct trios and apply MRPC
    
    
for (k in 1:id) {
      
# showing number of iterations every telliter
      
      
if( k %% telliter == 0 ) cat(paste("iteration", k, "complete\n"))
      
CNA.exp.meth <- matrix(0)
CNA <- NULL
exp <- NULL
meth <- NULL
      
CNA[k] <- paste0(full.range[k])
      
exp[k] <- paste0(full.range[k])
meth[k] <- paste0(full.range[k])

# Adding filters to skip when the data do not match
      
if(dim(n.n.data_CNA.2[which(n.n.data_CNA.2$Entrez_Gene_Id==CNA[k]),])[1]==0) next
if(dim(n.n.expression.2[which(n.n.expression.2$Entrez_Gene_Id==exp[k]),])[1]==0) next
      
# defining matrices to save results 
      

total.probes <- list()

adj.total <- matrix(0, dim(n.n.methylation.3[which(n.n.methylation.3$Entrez_Gene_Id==meth[k]),])[1], 9)
      
adj.cor <- matrix(0, dim(n.n.methylation.3[which(n.n.methylation.3$Entrez_Gene_Id==meth[k]),])[1], 16)
      
# Creating a loop to include repeated cases 
      
for (j in 1: dim(n.n.methylation.3[which(n.n.methylation.3$Entrez_Gene_Id==meth[k]),])[1]) {
        
        
if(dim(n.n.methylation.3[which(n.n.methylation.3$Entrez_Gene_Id==meth[k]),][j, ])[1]==0) next
        
        
CNA.exp.meth<-matrix(c(n.n.data_CNA.2[which(n.n.data_CNA.2$Entrez_Gene_Id==CNA[k]),],
                               n.n.expression.2[which(n.n.expression.2$Entrez_Gene_Id==exp[k]),],
                               n.n.methylation.3[which(n.n.methylation.3$Entrez_Gene_Id==meth[k]),][j, ]), dim(n.n.expression.2)[2] , 3)
        
      print(CNA.exp.meth[1:6, 1:3])
        
total.probes[j] <- CNA.exp.meth[1,3]
      

      
CNA.exp.meth[1:4, 1:3]
        
CNA.exp.meth.1 <- CNA.exp.meth[-1:-3,]
        
CNA.exp.meth.1[1:4, 1:3]  
      
CNA.exp.meth.2 <- matrix(as.numeric(unlist(CNA.exp.meth.1)), dim(CNA.exp.meth.1)[1], 3)
        
CNA.exp.meth.2[1:4, 1:3]
        

data <-CNA.exp.meth.2
        
        #print(data[1:4, 1:3])
        
dim(data)
        
n <- nrow (data)
        
        
colnames(data) <- c("V1", "T1","T2")
        
V <- colnames(data)     # Column names
        
#  print(cor(data, use = "complete.obs"))
        
cor12 <- cor(data[, 1], data[, 2], use = "complete.obs")
cor13 <- cor(data[, 1], data[, 3], use = "complete.obs")
cor23 <- cor(data[, 2], data[, 3], use = "complete.obs")

# Classical correlation
suffStat <- list(C = cor(data,use = "complete.obs"),
                         n = n)

# corrected correlation

suffStat$C[1,2] <- suffStat$C[2,1] <- cor12 
suffStat$C[1,3] <- suffStat$C[3,1] <- cor13 
suffStat$C[2,3] <- suffStat$C[3,2] <- cor23

cor.mtx <- suffStat$C
        
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
        
adj.total[j, ] <- matrix(Adj.infe1, 1, 9)
        

adj.cor[j, ] <- c(unlist(CNA.exp.meth[1, 3]), unlist(CNA.exp.meth[2, 3]),unlist(CNA.exp.meth[3, 3]), j,adj.total[j,], cor.mtx[1,2], cor.mtx[1,3], cor.mtx[2,3])
        
print(t(adj.cor[j, ]))
# write.table(t(adj.cor[j, ]), file="adj.cor.4.test.10.24.txt", sep = "\t", row.names = FALSE,col.names = FALSE , append =TRUE,quote=FALSE)
        
      }
    }
  }
}



sapply(full.range.2, function(x) f(x))


