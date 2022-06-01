top_ev <- function(data, info, data_type){
  
  #initialize the list
  mybiglist <- list()
  
  #loop through each column of PC
  for(i in 1:ncol(data)){
    
    #save one of the PCs to a variable
    pc <- data[,i]
  
    #ev <- sort(abs(pc), decreasing = TRUE)[1:10]
    
    #find the rows of the 10 highest absolute Eigenevector value for PC i
    row_num <- order(abs(pc), decreasing = TRUE)[1:10]
    
    #find the values in those rows
    ev <- data[row_num, i]
    
    #since the methylation data has an extra column with the Probe_ID
    #we will treat them separately
    if(data_type == "Methylation"){
    
    #find the gene name for the highest values
    gene_name <- info[row_num,ncol(info)]
    
    #find the probe_ID for the highest values
    probe_ID <- info[row_num,2]
    
    #merge the probe ID, gene name and highest value to res 
    res <- cbind(probe_ID, gene_name, ev)

    }else{
    
    #find the gene name for the highest values
    gene_name <- info[row_num,2]
    
    #merge the gene name and highest value to res 
    res <- cbind(gene_name, ev)
    
    }
    
    #save output from each loop to the list
    mybiglist[[i]] <- res

  }
 
  #return the list of top eigenvectors for all pc
  return(mybiglist)
  
}
