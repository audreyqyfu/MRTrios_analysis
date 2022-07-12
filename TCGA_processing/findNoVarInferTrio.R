no.var.gene <- function(gene.data){
  
  #find columns with no variance
  var.col <- apply(gene.data, 2, var, na.rm = TRUE)
  
  #return index of cols with 0 variance
  na.var.gene1 <- which(var.col == 0)
  
  #check if the column has NAs
  t <- apply(is.na(gene.data), 2, sum)
  
  #check if there are more 2 non-NA values
  na.var.gene2 <- which(t > (nrow(gene.data)-3))
  
  #combine the two missing columns info
  na.var.gene <- c(na.var.gene1, na.var.gene2)
  
  #return the results
  return(na.var.gene)
}


no.var.meth <- function(meth.data){
  
  #find columns with no variance
  var.col <- apply(meth.data, 2, var, na.rm = TRUE)
  
  #return index of cols with 0 variance
  na.var.meth1 <- which(var.col == 0)
  
  #check if the column has NAs
  t <- apply(is.na(meth.data), 2, sum)
  
  #check if there are more 2 non-NA values
  na.var.meth2 <- which(t > (nrow(meth.data)-3))
  
  #combine the two missing columns info
  na.var.meth <- c(na.var.meth1, na.var.meth2)
  
  #return the results
  return(na.var.meth)
  
}
