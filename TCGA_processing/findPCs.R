findPCs <- function(data, startCol, GeneNameCol, type.ind, com.ind, type){
  
  #find the rows in data that have atleast one NA values
  na.rows <- which(complete.cases(data[,startCol:ncol(data)]) == FALSE)
  
  #remove the the rows from the data
  data2 <- data[-na.rows,]
  dim(data2)
  data2[1:5,1:5]
  
  #finding common individuals between the 3 datasets and pos & neg ER individuals
  com.ind.type <- intersect(unlist(type.ind), com.ind)
  
  #find the column number of the individuals
  ind.col.data.type = match(com.ind.type, colnames(data2))
  
  #only save numeric values
  new.data.numeric <- t(data2[,ind.col.data.type])
  
  gene.names.data <- data2[,GeneNameCol]
  
  #get the columns with 0 variance and less than 3 non-NA values
  na.var.data2 <- no.var(new.data.numeric)
  
  if(length(na.var.data2) > 0){
    
    #remove those columns from the data and the gene names
    data.info <- gene.names.data[-na.var.data2]
    new.data.no.var <- new.data.numeric[,-na.var.data2]
  
  }else{
    
    #if there are no columns with no variance, keep the data as it is
    data.info <- gene.names.data
    new.data.no.var <- new.data.numeric
  }
  
  
  ############################################################################################
  #
  #                   Create the indices tables
  #
  ##########################################################################################
  
  # The first column is the rows numbers in the original data
  # The second column is the updated row numbers in the updated data
  # after removing NA values and no variance columns
  col1 <- 1:nrow(data)
  col2 <- rep(NA, nrow(data))
  to.rem <- c(na.rows, as.integer(names(na.var.data2)))
  col2[-to.rem] <- 1:(nrow(data) - length(to.rem))
  col.mtx.data <- cbind(col1,col2)
  
  
  ############################################################################################
  #
  #                   Get the PC scores
  #
  ##########################################################################################
  
  
  #calculate the PC score matrix
  pca.type <- prcomp(new.data.no.var, scale = TRUE)
  dim(pca.type$x) 
  
  
  #get the significant associated pcs
  data.with.conf = get.conf.matrix(new.data.no.var, pca.type$x, blocksize = 2000, apply.qval = TRUE)

 
  
  #replace the column numbers to corresponding gene names
  colnames(new.data.no.var) <- data.info
  dim(new.data.no.var)
  
  #return the 4 datasets
  return(list(pca.type$x, data.with.conf, col.mtx.data, new.data.no.var))
  
}
