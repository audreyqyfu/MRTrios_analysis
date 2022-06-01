library(data.table)
library(dplyr)

#Reading in the datasets
ev.neg.gene <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.eigenvectors.GeneExp.negER.new.txt"))
dim(ev.neg.gene)
ev.neg.gene[1:5,1:5]

ev.pos.gene <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.eigenvectors.GeneExp.posER.new.txt"))
dim(ev.pos.gene)
ev.pos.gene[1:5,1:5]

ev.neg.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.eigenvectors.Meth.negER.new.txt"))
dim(ev.neg.meth)
ev.neg.meth[1:5,1:5]

ev.pos.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/PCs/logit_data_PCs/PCA.eigenvectors.Meth.posER.new.txt"))
dim(ev.pos.meth)
ev.pos.meth[1:5,1:5]




#Gene Exp

#some pcs for Pos ER
top_ev(ev.pos.gene[,3:5], ev.pos.gene[,1:2], "Gene Exp")

#some pcs for Neg ER
top_ev(ev.neg.gene[,3:5], ev.neg.gene[,1:2], "Gene Exp")

#all pcs for Pos ER
top_ev(ev.pos.gene[,3:ncol(ev.pos.gene)], ev.pos.gene[,1:2], "Gene Exp")

#all pcs for Neg ER
top_ev(ev.neg.gene[,3:ncol(ev.neg.gene)], ev.neg.gene[,1:2], "Gene Exp")

#Meth

#some pcs for Pos ER
top_ev(ev.pos.meth[,3:5], ev.pos.meth[,c(1,2,ncol(ev.pos.meth))], "Methylation")

#some pcs for Neg ER
top_ev(ev.neg.meth[,3:5], ev.neg.meth[,c(1,2,ncol(ev.neg.meth))], "Methylation")

#all pcs for Pos ER
top_ev(ev.pos.meth[,3:(ncol(ev.pos.meth)-1)], ev.pos.meth[,c(1,2,ncol(ev.pos.meth))], "Methylation")

#all pcs for Neg ER
top_ev(ev.neg.meth[,3:(ncol(ev.neg.meth)-1)], ev.neg.meth[,c(1,2,ncol(ev.neg.meth))], "Methylation")
