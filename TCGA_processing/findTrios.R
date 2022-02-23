####################################
# function to find trios
# by matching methylation probes,
# gene expression, and copy number
# alteration (CNA)
####################################
# input
#
# methyl.data: methylation data matrix; probes in rows, individuals in columns
# expn.data: gene expression data matrix; genes in rows, individuals in columns
# cna.data: cna data matrix; genes in rows, individuals in columns
# entrez.list: master list of gene names with Entrez IDs
#
# output
# trios.rows: data matrix of 4 columns: gene name, row number
#             from each input dataset
#function <- findTrios (methyl.data, expn.data, cna.data, entrez.list) {
    # identify a list of unique genes
    # for methylation probes
    
    
    # identify trios by gene name (Hugo_Symbol)
    # call function trios()
    
    
    # identify additional trios by Entrez ID
    # call function entrez()
    
    
    # find additional matches using R package org.Hs.eg.db

    
    # return a data matrix for matched trios
    # in each row: gene name, row numbers from each dataset

#}




#function to match trios (works for one gene at a time)
trios <- function (gene.name) {
  
    #assign row numbers for the list data
    #so the genes in the same list with have the same row
    g1 <- rep(seq_along(meth.genes), sapply(meth.genes, length))
    g1[1:5]
    
    #match the gene name and fine the row
    meth.row <- g1[which(unlist(meth.genes) == gene.name)]
    meth.row
    
    
    ##### error here !!!!!!!!!!!!!!
    cna.row <- which(unlist(cna$Hugo_Symbol) == gene.name)
    cna.row

    if(length(cna.row) == 0){
      cna.row = NA
    }
    
    
    
    gene.row <- which(unlist(gene$Hugo_Symbol) == gene.name)
    gene.row
    
    if(length(gene.row) == 0){
      gene.row = NA
    }
    
    # gene.row <- match(gene.name, gene$Hugo_Symbol)
    # gene.row
    
        #we save the gene name, i (which is the row number in meth data)
        #and row numbers in cna & gene data
        results <- cbind(gene.name, meth.row, cna.row, gene.row)
        
        #return the result to the function
        return(results)
      
}




#function to find entrez ID rows (only works for the first match)
#entrez id are duplicates
entrez.cna <- function(cna.row){
  
  #use the cna row but take the column of the gene row
  #column 4 has gene rows in tmp and column 3 has cna rows in tmp
  entrez.row.cna <- which(cna$Entrez_Gene_Id == gene$Entrez_Gene_Id[as.integer(tmp[cna.row,4])])
  entrez.row.cna
  
  #check if the matched row is NA or not
  #If not proceed
  if(length(entrez.row.cna) > 0){
    
    #replace the matched rows with the respective NA spots in the tmp data
    tmp[cna.na.rows,3] <- entrez.row.cna
    
  }
  #return the data
  return(tmp)
 
}

entrez.gene <- function(gene.row){
  #use the gene row but take the column of the cna row
  entrez.row.gene <- which(gene$Entrez_Gene_Id) == cna$Entrez_Gene_Id[as.integer(tmp[gene.row,3])])
  entrez.row.gene
  
  #check if the matched row is NA or not
  #If not proceed
  if(length(entrez.row.gene) > 0){

    #replace the matched rows with the respective NA spots in the tmp data
    tmp[gene.na.rows,4] <- entrez.row.gene
  }
  
  #return the data
  return(tmp)
}
















