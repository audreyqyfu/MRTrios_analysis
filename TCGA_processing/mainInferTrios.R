#load the library
library(data.table)
library(splitstackshape)
library(MRGN, lib = "/mnt/ceph/kark6289/Rlibs")

#read in the datasets
gene.exp <- as.data.frame(fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt"))
dim(gene.exp)

cna <- as.data.frame(fread("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_CNA.txt"))
dim(cna)

TCGA.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt"))
dim(TCGA.meth)

#create a column with index
gene.exp$index <- 1:nrow(gene.exp)
dim(gene.exp)

TCGA.meth$index <- 1:nrow(TCGA.meth)
dim(TCGA.meth)

#counts the number of NAs in each row
#returns the rows that do not have a sum of NAs that match with the number of columns
#so basically removes the rows that have all NAs for methylation levels

gene <- gene.exp[rowSums(is.na(gene.exp[,3:(ncol(gene.exp)-1)])) != ncol(gene.exp[,3:(ncol(gene.exp)-1)]), ]
dim(gene)
gene[1:5,1:5]

meth <- TCGA.meth[rowSums(is.na(TCGA.meth[,5:(ncol(TCGA.meth)-1)])) != ncol(TCGA.meth[,5:(ncol(TCGA.meth)-1)]), ]
dim(meth)
meth[1:5,1:5]

#reading in the neg and pos ER individuals data
neg.ind <- fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.neg.patient.txt", header = FALSE)
dim(neg.ind)

pos.ind <- fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.pos.patient.txt", header = FALSE)
dim(pos.ind)

#finding common individuals among the 3 datasets
match1 = intersect(colnames(gene)[3:ncol(gene)], colnames(meth)[5:ncol(meth)]) #length 787
com.ind = intersect(colnames(cna)[3:ncol(cna)], match1) #length 777

#reading in the Trios data
trios <- data.frame(fread("/mnt/ceph/kark6289/test_trio/trios/trio.final.protein.coding.txt"))

#finding common individuals between the 3 datasets and pos & neg ER individuals
com.ind.pos <- intersect(unlist(pos.ind[,2]), com.ind)
com.ind.neg <- intersect(unlist(neg.ind[,2]), com.ind)

#reading in the PC score matrix
pc.gene.pos = as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.GeneExp.posER.txt"))
pc.meth.pos = as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.Meth.posER.txt"))

pc.gene.neg = as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.GeneExp.negER.txt"))
pc.meth.neg = as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.Meth.negER.txt"))

#reading in the sig asso pcs
gene.sig.asso.pcs.pos <- readRDS("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/gene.posER.sig.asso.pcs.RData")
meth.sig.asso.pcs.pos <- readRDS("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/meth.posER.sig.asso.pcs.RData")

gene.sig.asso.pcs.neg <- readRDS("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/gene.negER.sig.asso.pcs.RData")
meth.sig.asso.pcs.neg <- readRDS("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/meth.negER.sig.asso.pcs.RData")


mrgn.InferTrio <- function(trios, gene.exp, TCGA.meth, gene, meth, pc.gene, pc.meth, com.ind.type, gene.sig.asso.pcs, meth.sig.asso.pcs, type.ER){
  
  #only save numeric values
  new.gene.nona <- t(gene[,-c(1,2,ncol(gene))])
  new.meth.nona <- t(meth[,-c(1,2,3,4,ncol(meth))])
  

  #match the individuals between data and pc
  ind.com.pc.trio.gene <- match(unlist(pc.gene[,1]), rownames(new.gene.nona))
  ind.com.pc.trio.meth <- match(unlist(pc.meth[,1]), rownames(new.meth.nona))

  #use the common individuals
  gene.data <- new.gene.nona[ind.com.pc.trio.gene,]
  meth.data <- new.meth.nona[ind.com.pc.trio.meth,]

  na.var.gene <- no.var.gene(gene.data)
  na.var.meth <- no.var.gene(meth.data)

  if(length(na.var.gene)>0){
    
  #create a matrix with column index of the data
  #with new column index in another column after removal of columns with 0 variance
  col1 <- 1:ncol(gene.data)
  col2 <- rep(NA, ncol(gene.data))
  col2[-na.var.gene] <- 1: (ncol(gene.data) - length (na.var.gene))
  col.mtx.gene <- cbind (col1, col2)
  }
  
  if(length(na.var.meth)>0){
    
  #create a matrix with column index of the data
  #with new column index in another column after removal of columns with 0 variance
  col1 <- 1:ncol(meth.data)
  col2 <- rep(NA, ncol(meth.data))
  col2[-na.var.meth] <- 1: (ncol(meth.data) - length (na.var.meth))
  col.mtx.meth <- cbind (col1, col2)
  }
  
   
 
    #matching individuals between common positive individuals and respective datasets
    ind.col.meth = match(com.ind.type, colnames(TCGA.meth))
    ind.col.gene = match(com.ind.type, colnames(gene.exp))
    ind.col.cna = match(com.ind.type, colnames(cna))
    
    #match individuals
    com.ind.trio.pc.gene <- match(colnames(gene.exp[,ind.col.gene]), pc.gene[,1])
    com.ind.trio.pc.meth <- match(colnames(TCGA.meth[,ind.col.meth]), pc.meth[,1])

    result = NULL
    
    for(i in 1:nrow(trios)){
    
    if(is.na(trios[i,3]) == FALSE & is.na(trios[i,4]) == FALSE) {
      
      if(rowSums(is.na(gene.exp[trios[i,4],3:(ncol(gene.exp)-1)])) != ncol(gene.exp[,3:(ncol(gene.exp)-1)])){
      
        
        #create the trio
        trio.cna = t(cna[trios[i,3],ind.col.cna])
        trio.gene = t(gene.exp[trios[i,4],ind.col.gene])
        trio.meth = t(TCGA.meth[trios[i,2],ind.col.meth])
        
        trio.mat = cbind(trio.cna, trio.gene, trio.meth)
        
        #finding the gene row after removal of NAs
        gene.row.nona <- which(gene$index == gene.exp$index[as.numeric(colnames(trio.mat)[2])])
        meth.row.nona <- which(meth$index == TCGA.meth$index[as.numeric(colnames(trio.mat)[3])])
        
        if(length(na.var.gene) > 0){
          
          gene.row.novar <- col.mtx.gene[gene.row.nona,2]
          
        }else{
          gene.row.novar <- gene.row.nona
        }
        
        if(length(na.var.meth) > 0){
          
          meth.row.novar <- col.mtx.meth[meth.row.nona,2]
          
        }else{
          
          meth.row.novar <- meth.row.nona
          
        }
        
        #find common pcs between gene exp and meth
        #com.sig.asso.pcs <- union(unlist(gene.sig.asso.pcs[gene.row.nona,]), unlist(meth.sig.asso.pcs[meth.row.nona,]))
        
        #get the sig pcs from the pc score matrix
        sig.pc.gene <- pc.gene[com.ind.trio.pc.gene,(unlist(gene.sig.asso.pcs[gene.row.novar,])+1)]
        sig.pc.meth <- pc.meth[com.ind.trio.pc.meth,(unlist(meth.sig.asso.pcs[meth.row.novar,])+1)]
        
        total.pc.count <- length(unlist(gene.sig.asso.pcs[gene.row.novar,])) + length(unlist(meth.sig.asso.pcs[meth.row.novar,]))
        
        #create matrix
        final.mat <- cbind(trio.mat, sig.pc.gene, sig.pc.meth)
        
        #infer the trio
        res = infer.trio(as.data.frame(final.mat), use.perm = FALSE)
        
        final <- cbind(rownames(trios[i,]), res$Inferred.Model, total.pc.count)
        
      }
    }
      
      #save the result from each loop
      result <- rbind(result,final)
      
    }
    
        print(type.ER)
    
        colnames(result) <- c("index", "Model_type", "Total_PC_count")
        return(result)
    
    
        #write to a file
        #write.table(final, file = "/mnt/ceph/kark6289/TCGA_analysis/MRGN_InferTrio/model.trio.MRGN.txt", sep = "\t", row.names = FALSE,
         #           col.names = FALSE, append = TRUE, quote=FALSE)

}

#If extracting the results for Positive ER type
mrgn.InferTrio(trios[1:3,], gene.exp, TCGA.meth, gene, meth, pc.gene.pos, pc.meth.pos, com.ind.pos, gene.sig.asso.pcs.pos, meth.sig.asso.pcs.pos, type.ER = "Positive")

#If extracting the results for Negative ER type
mrgn.InferTrio(trios[1:3,], gene.exp, TCGA.meth, gene, meth, pc.gene.neg, pc.meth.neg, com.ind.neg, gene.sig.asso.pcs.neg, meth.sig.asso.pcs.neg, type.ER = "Negative")

