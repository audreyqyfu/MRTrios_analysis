#function to use MRGN
analyzeTrios <- function(TCGA.meth, gene.exp, cna, trios, pc.meth, pc.gene, meth.sig.asso.pcs, gene.sig.asso.pcs, clinical, meth.table, gene.table, path){
  
  # find the common individuals between the 3 datasets
  # pc matrix has common individuals from meth and gene exp
  com.ind <- intersect(rownames(pc.meth),colnames(cna))
  
  #find the rows in clinical data for the common individuals
  rows.clinical <- match(com.ind, unlist(clinical[,1]))
  
  #extract the age and race for those individuals
  age <- clinical.pos[rows.clinical,2]
  race <- clinical.pos[rows.clinical,3]
  
  #find the rows for the common individuals in the resp datasets
  ind.col.cna <- match(com.ind, colnames(cna))
  ind.col.gene <- match(com.ind, colnames(gene.exp))
  ind.col.meth <- match(com.ind, colnames(TCGA.meth))
  
  #initialize tmp
  tmp <- NULL
  
  #begin the loop for rows in trios
  for(i in 270000:nrow(trios)){
    
    #check if the values for cna and gene exp are NA or not
    if(is.na(trios[i,3]) == FALSE & is.na(trios[i,4]) == FALSE) {
      
      if(rowSums(is.na(gene.exp[trios[i,4],3:ncol(gene.exp)])) != ncol(gene.exp[,3:ncol(gene.exp)])){
        
        #print the row number
        print(i)
        
        #extract data for each dataset and create the trio
        trio.cna = t(cna[trios[i,3],ind.col.cna])
        trio.gene = t(gene.exp[trios[i,4],ind.col.gene])
        trio.meth = t(TCGA.meth[trios[i,2],ind.col.meth])
        
        #combine the matrix
        trio.mat = cbind(trio.cna, trio.gene, trio.meth)
        
        #extract the updated row number from the indices table
        row.sig.pcs.meth <- unlist(meth.table[trios[i,2],2])
        row.sig.pcs.gene <- unlist(gene.table[trios[i,4],2])
        
        #find the row numbers for the common individuals in the pc score matrix
        com.ind.meth <- match(com.ind, rownames(pc.meth))
        com.ind.gene <- match(com.ind, rownames(pc.meth))
        
        #extract the column numbers for sig pcs
        sig.pcs.cols.meth <- as.integer(unlist(meth.sig.asso.pcs[[1]][[row.sig.pcs.meth]]))
        sig.pcs.cols.gene <- as.integer(unlist(gene.sig.asso.pcs[[1]][[row.sig.pcs.gene]]))
        
        #remvove PCs greater than 50
        insig.meth <- which(sig.pcs.cols.meth > 50)
        
        if(length(insig.meth) > 0){
          
          sig.pcs.cols.meth <- sig.pcs.cols.meth[-insig.meth]
          
        }else{
          
          sig.pcs.cols.meth <- sig.pcs.cols.meth
          
        }
        
        
        insig.gene <- which(sig.pcs.cols.gene > 50)
        
        if(length(insig.gene) > 0){
          
          sig.pcs.cols.gene <- sig.pcs.cols.gene[-insig.gene]
          
        }else{
          
          sig.pcs.cols.gene <- sig.pcs.cols.gene
          
        }
        
        #extract the sig columns from the pc matrix with the common individuals  
        sig.pc.gene <- pc.gene[com.ind.gene, sig.pcs.cols.gene]
        sig.pc.meth <- pc.meth[com.ind.meth, sig.pcs.cols.meth]
        
        #count the total number of sig pcs
        total.pc.count <- length(unlist(gene.sig.asso.pcs[[1]][[row.sig.pcs.gene]])) + length(unlist(meth.sig.asso.pcs[[1]][[row.sig.pcs.meth]]))
        
        #create matrix with the trios and the confounding variables
        final.mat <- cbind(trio.mat, sig.pc.gene, sig.pc.meth, age, race)
        
        #apply MRGN and infer the trio
        res = infer.trio(as.data.frame(final.mat), use.perm = TRUE, is.CNA = TRUE, nperms = 500)
        
        #combine the row number of trios, model type, and pc count
        final <- cbind(i, res, total.pc.count)
        
        #combine the value for each trio in rows
        tmp <- rbind(tmp, final)
        
        #write to a file
        write.table(final, file = path, sep = "\t", row.names = FALSE,
                     col.names = FALSE, append = TRUE, quote=FALSE)
      
      } 
    }
  }
  
  #return the dataset
  return(tmp)
  
}
