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


#function to match trios (works for one gene at a time)
trios <- function (gene.name, meth, cna, gene, meth.genes, dups.cna, dups.gene, na.meth, na.gene) {

    #assign row numbers for the list data
    #so the genes in the same list with have the same row
    g1 <- rep(seq_along(meth.genes), sapply(meth.genes, length))
    g1[1:5]

    #match the gene name and find the row in methylation data
    meth.row <- g1[which(unlist(meth.genes) == gene.name)]
    meth.row

    #match the gene name and find the row in CNA data
    cna.row <- which(unlist(cna$Hugo_Symbol) == gene.name)
    cna.row

    #if no match found assign NA
    if(length(cna.row) == 0){
      cna.row = NA
    }

    #match the gene name and find the row in Gene Exp data
    gene.row <- which(unlist(gene$Hugo_Symbol) == gene.name)
    gene.row

    #if no match found assign NA
    if(length(gene.row) == 0){
      gene.row = NA
    }

    #initialize results
    results = NULL

        #since there are multiple matches for each gene name in the meth data
        #we loop through each row and make sure it does not have NA for every individual
        for (i in 1:length(meth.row)){
        #print(i)

        #check for NA rows
        if((as.integer(meth.row[i]) %in% as.integer(na.meth)) == FALSE & (as.integer(gene.row) %in% as.integer(na.gene)) == FALSE) {

    #print(i)
 #we save the gene name, i (which is the row number in meth data)
    #and row numbers in cna & gene data
    results <- rbind(results,cbind(gene.name, meth.row[i], cna.row, gene.row))
}
}
    #return the result to the function
    return(results)
}



#function to find entrez ID rows in CNA data
entrez.cna <- function(cna.row, data, cna, gene){

  #use the cna row but take the column of the gene row
  #column 4 has gene rows in tmp and column 3 has cna rows in tmp
  entrez.row.cna <- which(cna$Entrez_Gene_Id == gene$Entrez_Gene_Id[as.integer(data[cna.row,4])])
  entrez.row.cna

  #check if the matched row is NA or not
  #If not proceed
  if(length(entrez.row.cna) < 1){

    #replace the matched rows with the respective NA spots in the tmp data
    entrez.row.cna <- NA
  }
  #return the data
  return(entrez.row.cna)

}

#function to find entrez ID rows in Gene Exp data
entrez.gene <- function(gene.row, data, cna, gene){

   #use the gene row but take the column of the cna row
  entrez.row.gene <- which(gene$Entrez_Gene_Id == cna$Entrez_Gene_Id[as.integer(data[gene.row,3])])
  entrez.row.gene

  #check if the matched row is NA or not
  #If not proceed
  if(length(entrez.row.gene) < 1){

    #replace the matched rows with the respective NA spots in the tmp data
    entrez.row.gene <- NA
  }
  #return the data
  return(entrez.row.gene)

}
