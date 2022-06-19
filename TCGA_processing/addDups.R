cna.add.dups <- function(final, cna){
 
  new.row <- NULL

  #rows with duplicates in CNA data
dup.CNA <- findDupsCNA(cna)
  
 #which rows have NA in the cna column in trio file
na.row.cna <- which(is.na(final[,3]) == TRUE)

#load the r package
xx <- as.list(org.Hs.egALIAS2EG)

#find the unique genes that have NA in cna column
genes.na.row <- unique(unlist(final[na.row.cna,1]))
dup.final <- NULL


for(i in 1:length(genes.na.row)){
  print(i)

  #find the gene name for the missing entrez id
  gene.match <- which(names(xx) == as.character(genes.na.row[i]))
  gene.match

  if(length(gene.match) > 0){

    #find the entrez id
    entz.id.new <- xx[[gene.match]]
    entz.id.new

    if(length(entz.id.new) == 1){

      #find the row number for that entrez id in cna data
      new.row <- which(cna$Entrez_Gene_Id == entz.id.new)
      new.row

    }else if(length(entz.id.new) > 1){
      new.row = NULL

      for(j in 1:length(entz.id.new)){

        row <- which(cna$Entrez_Gene_Id == entz.id.new[j])
        new.row <- c(new.row, row)

      }
    }

     if(length(new.row) > 0){
      
      save.row <- NULL
      
      for(i in 1:length(new.row)){
        
        if((new.row[i] %in% dup.CNA) == FALSE){
          save.row = c(save.row,new.row[i])
        }
      }
    
    new.row <- save.row
   }
     
 #check if length greater than 0
    if(length(new.row) == 1){
      print(i)

        #find which rows have the gene in cna data
        final.cna.row <- which(final[,1] == as.character(genes.na.row[i]))
        final[final.cna.row,3] <- new.row
    }

    #check if length greater than 1
    if(length(new.row) > 1){
      print(i)


      #find the rows in the cna data that contain the specified gene
      final.cna.row <- which(final[,1] == as.character(genes.na.row[i]))

      # we save all the remaining values (except the first one) to dup.final
        for(k in 2:length(new.row)){

          dup.final <- rbind(dup.final,final[final.cna.row,])
          dup.final[which(is.na(dup.final[,3]) == TRUE),3] <- new.row[k]

        }

        # and save the first value to the original file
        final[final.cna.row,3] <- new.row[1]

    }

  }
}

return(list(final,dup.final))

}


gene.add.dups <- function(final, dup.final, gene){

  new.row <- NULL

  #rows with duplicates in Gene Exp data
dup.GENE <- findDupsGENE(gene)
  
  #which rows have NA in the gene exp column in trio file
  na.row.gene <- which(is.na(final[,4]) == TRUE)

  #load the r package
  xx <- as.list(org.Hs.egALIAS2EG)

  #find the unique genes that have NA in gene exp column
  genes.na.row <- unique(unlist(final[na.row.gene,1]))

for(i in 1:length(genes.na.row)){
  print(i)

  #find the gene name for the missing entrez id
  gene.match <- which(names(xx) == as.character(genes.na.row[i]))
  gene.match

  #check if there is match in the R library for the specified gene name
  if(length(gene.match) > 0){

  #find the entrez id
  entz.id.new <- xx[[gene.match]]
  entz.id.new

  #check the length of the entz id and find the corresponding row in the data
  if(length(entz.id.new) == 1){

    #find the row number for that entrez id in cna data
    new.row <- which(gene$Entrez_Gene_Id == entz.id.new)
    new.row

   #check if length is greater than 1
  }else if(length(entz.id.new) > 1){
    new.row = NULL

    #loop through the entz id to find the match in gene exp data
    for(j in 1:length(entz.id.new)){

      #find the match for each entz id
      row <- which(gene$Entrez_Gene_Id == entz.id.new[j])

      #save every row to new.row
      new.row <- c(new.row, row)

    }
  }
 
   if(length(new.row) > 0){
      
      save.row <- NULL
      
    for(i in 1:length(new.row)){
      
        if((new.row[i] %in% dup.GENE) == FALSE){
          save.row = c(save.row,new.row[i])
        }
      }
    
    new.row <- save.row
   }
   
  #check if length greater than 0
  if(length(new.row) == 1){


    #find which rows have the gene in gene.exp data
    final.gene.row <- which(final[,1] == as.character(genes.na.row[i]))
    final[final.gene.row,4] <- new.row

  }

  #check if length greater than 1
  if(length(new.row) > 1){

    #find the rows in the gene exp data that contain the specified gene
    final.gene.row <- which(final[,1] == as.character(genes.na.row[i]))

    #finding the rows for the gene in dup.final since we already created this data
    # while finding rows for cna data
    dup.gene.row <- which(dup.final[,1] == as.character(genes.na.row[i]))

    #length decides if the gene is already in dup.final or not
    # if the length is greater than 0, we replace NA for gene.row in the existing data
    if(length(dup.gene.row) > 0){

      # save the first value to the original file
      final[final.gene.row,4] <- new.row[1]

      #for each new.row value starting from 2 save to dup.final
      for(k in 2:length(new.row)){

      dup.final[dup.gene.row[1:length(final.gene.row)],4] <- new.row[k]

      }
  # if it is not we create dup.final for the gene
    }else{

      #for each new.row value starting from 2 save to dup.final
      for(k in 2:length(new.row)){

      #copy the corresponding data from final to dup.final
      dup.final <- rbind(dup.final,final[final.gene.row,])

      #replace the NA with the new.row in gene.row column
      dup.final[which(is.na(dup.final[,4]) == TRUE),4] <- new.row[k]

      }

      # save the first value to the original file
      final[final.gene.row,4] <- new.row[1]

    }
   }
  }
}

  #combine the datasets into one
  final.res = rbind(final,dup.final)

  #rename the second column
  colnames(final.res)[2] = "meth.row"

  return(final.res)


}
