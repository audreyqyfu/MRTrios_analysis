findTrioAll <- function(meth.data, cna.data, gene.data) {

#rows with duplicates in CNA data
dup.CNA <- findDupsCNA(cna.data)

#rows with duplicates in Gene Exp data
dup.GENE <- findDupsGENE(gene.data)

#rows with NAs for all individual in Methylation data
na.METH <- removeNA.meth(meth.data)

#rows with NAs for all individual in Gene Exp data
na.GENE <- removeNA.gene(gene.data)

#since the data has multiple genes in one row separated by ";", we split them
split.genes <- unlist(strsplit(as.character(meth.data$Gene_Symbol), ';'))

#then find the unique genes from the methylation data
uni <- na.omit(unique(split.genes))

#remove any empty values in the uni list and skip the duplicated rows in CNA data
uni2 <- uni[-na.omit(match(cna$Hugo_Symbol[dup.CNA], uni))]

#skip the duplicated rows in Gene Exp data
unique.genes <- uni2[-na.omit(match(gene$Hugo_Symbol[dup.GENE], uni2))]

#to find the row in the methylation data
meth.genes <- strsplit(as.character(meth.data$Gene_Symbol), ';')

#initialize the tmp variable
tmp <- NULL

    #loop through each gene
    for(i in 1:length(unique.genes)){
    print(i)

        #apply the function with the provided data
        tmp <- rbind(tmp,trios(unique.genes[i], meth.data, cna.data, gene.data, meth.genes, dup.CNA, dup.GENE, na.METH, na.GENE))

  }

# dim(tmp)
#print(tmp[1:5,])

#find the rows that have NAs for each data type
meth.na.rows <- which(is.na(tmp[,2]))
cna.na.rows <- which(is.na(tmp[,3]))
gene.na.rows <- which(is.na(tmp[,4]))

#loop through the rows that did not have a gene match in CNA data
for(i in 1:length(cna.na.rows)){

    #find the entrez id match and replace the NA
    tmp[cna.na.rows[i],3] <- entrez.cna(cna.na.rows[i], tmp, cna.data, gene.data)
}

print(tmp[1:5,])

#loop through the rows that did not have a gene match in Gene Exp data
for(i in 1:length(gene.na.rows)){

    #find the entrez id match and replace the NA
    tmp[gene.na.rows[i],4] <- entrez.gene(gene.na.rows[i], tmp, cna.data, gene.data)
}

print(tmp[1:5,])

#return the final data
return(tmp)

}
