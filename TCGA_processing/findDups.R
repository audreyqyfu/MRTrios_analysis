findDups <- function(data){

  # Finding duplicates in entrez id

  #gives rows for which there are no NA as entrez id
  nona.ent <- which(!is.na(data$Entrez_Gene_Id))

  #the entrez ids of the above rows
  nona.entid <- data$Entrez_Gene_Id[nona.ent]
  nona.entid[1:5]

  #checks if the entred ids have duplicates or not (TRUE/FALSE)
  d.ent <- duplicated(nona.entid)

  #find which rows have duplicates (gives only the first one not the appearances)
  dups.rows <- which(d.ent == TRUE)

  #pick the first duplicated row and find the entrez id
  entz.ids <- nona.entid[dups.rows]

  #initialize variable for duplicated rows in CNA data
  dups.row.data <- NULL

  #find the rows where this particular entrez id is duplicated
  for(i in 1:length(entz.ids)){

  #since there are multiple entrez ids, we loop through each one of them
  dups <- which(data$Entrez_Gene_Id == entz.ids[i])

  #save all the row numbers
  dups.row.data <- c(dups.row.data, dups)
  }

  #return the result
  return(dups.row.data)
}




findDupsCNA <- function(cna){

  # Finding duplicates in entrez id

  #gives rows for which there are no NA as entrez id
  nona.ent <- which(!is.na(cna$Entrez_Gene_Id))

  #the entrez ids of the above rows
  nona.entid <- cna$Entrez_Gene_Id[nona.ent]
  nona.entid[1:5]

  #checks if the entred ids have duplicates or not (TRUE/FALSE)
  d.ent <- duplicated(nona.entid)

  #find which rows have duplicates (gives only the first one not the appearances)
  dups.rows <- which(d.ent == TRUE)

  #pick the first duplicated row and find the entrez id
  entz.ids <- nona.entid[dups.rows]

  #initialize variable for duplicated rows in CNA data
  dups.row.cna <- NULL

  #find the rows where this particular entrez id is duplicated
  for(i in 1:length(entz.ids)){

  #since there are multiple entrez ids, we loop through each one of them
  dups <- which(cna$Entrez_Gene_Id == entz.ids[i])

  #save all the row numbers
  dups.row.cna <- c(dups.row.cna, dups)
  }

  #return the result
  return(dups.row.cna)
}


findDupsGENE <- function(gene){

  #gives rows for which there are no NA as entrez id
  nona.ent <- which(!is.na(gene$Entrez_Gene_Id))

  #the entrez ids of the above rows
  nona.entid <- gene$Entrez_Gene_Id[nona.ent]
  nona.entid[1:5]

  #checks if the entred ids have duplicates or not (TRUE/FALSE)
  d.ent <- duplicated(nona.entid)

  #find which rows have duplicates (gives only the first one not the appearances)
  dups.rows <- which(d.ent == TRUE)

  #pick the first duplicated row and find the entrez id
  entz.ids <- nona.entid[dups.rows]

  #initialize variable for duplicated rows in
  dups.row.gene <- NULL

  #find the rows where this particular entrez id is duplicated
  for(i in 1:length(entz.ids)){

  #since there are multiple entrez ids, we loop through each one of them
  dups <- which(gene$Entrez_Gene_Id == entz.ids[i])

  #save all the row numbers
  dups.row.gene <- c(dups.row.gene, dups)
  }

   #return the result
   return(dups.row.gene)
}
