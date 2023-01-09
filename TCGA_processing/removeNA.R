removeNA <- function(data){

   #find the rows that have NA as values for all the individuals in Methylation data
   na.cols.meth <- which(rowSums(is.na(TCGA.meth[,nStart:ncol(TCGA.meth)])) == ncol(TCGA.meth[,nStart:ncol(TCGA.meth)]))

   #return the result
   return(na.cols)

}


removeNA.meth <- function(TCGA.meth){

   #find the rows that have NA as values for all the individuals in Methylation data
   na.cols.meth <- which(rowSums(is.na(TCGA.meth[,5:ncol(TCGA.meth)])) == ncol(TCGA.meth[,5:ncol(TCGA.meth)]))

   #return the result
   return(na.cols.meth)


}


removeNA.gene <- function(gene.exp){

   #find the rows that have NA as values for all the individuals in Gene Exp data
   na.cols.gene <- which(rowSums(is.na(gene.exp[,3:ncol(gene.exp)])) == ncol(gene.exp[,3:ncol(gene.exp)]))

   #return the result
   return(na.cols.gene)


}
