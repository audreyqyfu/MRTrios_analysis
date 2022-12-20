combined_df1 <- read.delim("/mnt/ceph/kark6289/PCandTrioAnalysis/output.09.30/model.trio.MRGN.all.posER", header = FALSE)

combined_df2 <- read.delim("/mnt/ceph/kark6289/PCandTrioAnalysis/output.09.30/model.trio.MRGN.all.negER", header = FALSE)

res <- combined_df1  
#res <- combined_df2

#remove rows that have NA
res <- res[-which(is.na(res), arr.ind=TRUE),]

#rename the column names
colnames(res) <- c("trio_row", "b11", "b12", "b21", "b22", "V1:T1", "V1:T2", "pb11", "pb12", "pb21", "pb22", "pV1:T1", "pV1:T2", "Minor.freq", "Inferred.Model", "PC_count")

#create a new column for new model name after filtering
res$Inferred.Model2[1:nrow(res)] <- rep(NA)

#loop through all the rows in the data
for(i in 1:nrow(res)){

  print(i)
  
#if classification is either M0, M1, or M2
if(res$Inferred.Model[i] == "M0.1" | res$Inferred.Model[i] == "M0.2" | res$Inferred.Model[i] == "M2.1" | res$Inferred.Model[i] == "M2.2"){
  
  #we check the marginal p values and if they are less than 0.01, we re classify as "Other"
  if(res$`pV1:T2`[i] < 0.01 & res$`pV1:T1`[i] < 0.01){
    
    res$Inferred.Model2[i] <- "Other"
    
  }
 }else{
  
  res$Inferred.Model2[i] <- res$Inferred.Model[i]
  
 }
  
  #for all the models, we check the marginal pvalues, if they are greather than 0.05, we re classify as "Other"
  if(res$`pV1:T2`[i] > 0.05 & res$`pV1:T1`[i] > 0.05){
    
    res$Inferred.Model2[i] <- "Other"
    
  }else{
    
    res$Inferred.Model2[i] <- res$Inferred.Model[i]
    
  }
  
}


# 292455     17
# 292331     17

#write.table(res, file = "/mnt/ceph/kark6289/PCandTrioAnalysis/output.09.30/model.trio.MRGN.all.posER.reclassify.txt", sep = "\t", row.names = TRUE,
#            col.names = NA, append = TRUE, quote=FALSE)

#write.table(res, file = "/mnt/ceph/kark6289/PCandTrioAnalysis/output.09.30/model.trio.MRGN.all.negER.reclassify.txt", sep = "\t", row.names = TRUE,
#            col.names = NA, append = TRUE, quote=FALSE)
