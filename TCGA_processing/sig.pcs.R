####################################################

           THIS ONLY HAS CODE FOR METHYLATION NEGATIVE ER FOR NOW
           
##########################################################



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

#counts the number of NAs in each row
#returns the rows that do not have a sum of NAs that match with the number of columns
#so basically removes the rows that have all NAs for methylation levels

gene <- gene.exp[rowSums(is.na(gene.exp[,3:ncol(gene.exp)])) != ncol(gene.exp[,3:ncol(gene.exp)]), ]
dim(gene)
gene[1:5,1:5]

meth <- TCGA.meth[rowSums(is.na(TCGA.meth[,5:899])) != ncol(TCGA.meth[,5:899]), ]
dim(meth)
meth[1:5,1:5]

#read in the pc score matrix
  pc.meth = as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.Meth.negER.txt"))
  
  #only save numeric values
  new.meth.nona <- t(meth[,-(1:4)])
  
  #match the individuals between data and pc
  ind.com.pc.trio <- match(unlist(pc.meth[,1]), rownames(new.meth.nona))
  
  #use the common individuals
  meth.data <- new.meth.nona[ind.com.pc.trio,]
  
 #find columns with no variance
  var.col <- apply(meth.data, 2, var, na.rm = TRUE)
  na.var <- which(var.col == 0)
  
  #check length of na.var
  if(length(na.var) > 0){
    
  meth.with.conf = get.conf(meth.data[,-unname(na.var)], pc.meth[,-1], blocksize = 2000, return.for.trios = FALSE, method = "correlation")
  
  }else{
    
    meth.with.conf = get.conf(meth.data, pc.meth[,-1], blocksize = 2000, return.for.trios = FALSE, method = "correlation")
    
  }
  
  #breakdown the nested list into a single list
  tmp <- as.data.frame(do.call(rbind, meth.with.conf$sig.asso.pcs))
  
  #first set of sig pcs
  tmp[1,]
  
  #save the data
  #saveRDS(tmp, file="/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/meth.negER.sig.asso.pcs.RData")
  
  #read in the save data
  #sig.pc <- readRDS("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/meth.negER.sig.asso.pcs.RData")
  
  
