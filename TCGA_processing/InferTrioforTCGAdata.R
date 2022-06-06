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

#counts the number of NAs in each row
#returns the rows that do not have a sum of NAs that match with the number of columns
#so basically removes the rows that have all NAs for methylation levels

gene <- gene.exp[rowSums(is.na(gene.exp[,3:ncol(gene.exp)])) != ncol(gene.exp[,3:ncol(gene.exp)]), ]
dim(gene)
gene[1:5,1:5]

meth <- TCGA.meth[rowSums(is.na(TCGA.meth[,5:899])) != ncol(TCGA.meth[,5:899]), ]
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
trios <- data.frame(fread("/mnt/ceph/kark6289/test_trio/trios/Trios.final2.txt"))

#finding common individuals between the 3 datasets and pos & neg ER individuals
com.ind.pos <- intersect(unlist(pos.ind[,2]), com.ind)
com.ind.neg <- intersect(unlist(neg.ind[,2]), com.ind)

#reading in the PC score matrix
  pc.gene = as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.GeneExp.posER.txt"))
  pc.meth = as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.Meth.posER.txt"))
  
  #reading in the sig asso pcs
  gene.sig.asso.pcs <- readRDS("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/gene.posER.sig.asso.pcs.RData")
  meth.sig.asso.pcs <- readRDS("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/meth.posER.sig.asso.pcs.RData")
  
  #matching individuals between common positive individuals and respective datasets
  ind.col.meth = match(com.ind.pos, colnames(meth))
  ind.col.gene = match(com.ind.pos, colnames(gene))
  ind.col.cna = match(com.ind.pos, colnames(cna))
  
  #mybiglist <- list()
  
  for(i in 1:200){
    
    #create the trio
    trio.cna = t(cna[trios[i,3],ind.col.cna])
    trio.gene = t(gene.exp[trios[i,4],ind.col.gene])
    trio.meth = t(TCGA.meth[trios[i,2],ind.col.meth])
    
    trio.mat = cbind(trio.cna, trio.gene, trio.meth)
    
    #match individuals
    com.ind.trio.pc.gene <- match(rownames(trio.mat), pc.gene[,1])
    
    com.ind.trio.pc.meth <- match(rownames(trio.mat), pc.meth[,1])
    
    #finding the gene row after removal of NAs
    gene.row.nona <- which(gene$Hugo_Symbol == gene.exp$Hugo_Symbol[as.numeric(colnames(trio.mat)[2])])
    meth.row.nona <- which(meth$Row.names == TCGA.meth$Row.names[as.numeric(colnames(trio.mat)[3])])
    
    
    #find common pcs between gene exp and meth
    com.sig.asso.pcs <- union(unlist(gene.sig.asso.pcs[gene.row.nona,]), unlist(meth.sig.asso.pcs[meth.row.nona,]))
    
    #get the sig pcs from the pc score matrix
    sig.pc.gene <- pc.gene[com.ind.trio.pc.gene,(com.sig.asso.pcs+1)]
    
    sig.pc.meth <- pc.meth[com.ind.trio.pc.meth,(com.sig.asso.pcs+1)]
    
    #create matrix
    final.mat <- cbind(trio.mat, sig.pc.gene, sig.pc.meth)
    
  }
  
  
