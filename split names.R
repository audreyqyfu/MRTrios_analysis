#set the directory
setwd("/mnt/ceph/megheib/Cancer_Genomics/TCGA_data_results")
getwd()

#load the library
library(data.table)
library(splitstackshape)

load("TCGA.meth.RData")
dim(TCGA.meth)
meth = TCGA.meth

for(i in 5:ncol(meth)){

  print(i)

  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(colnames(meth)[i]), "[.]"))
  split.ID

  #now merge the first 4 parts
  new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3],"-",split.ID[4], sep = "")
  new.rowname

  #replace the old column name to the new column name
  colnames(meth)[i] <- new.rowname


}

#save it to a new txt file
write.table(meth, file = "/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.txt", sep = "\t", row.names = FALSE,
            col.names = TRUE, append = TRUE, quote=FALSE)


library(data.table)
TCGA.meth <- as.data.frame(fread("/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt"))
dim(TCGA.meth)

data <- TCGA.meth[rowSums(is.na(TCGA.meth[,5:ncol(TCGA.meth)])) != ncol(TCGA.meth[,5:ncol(TCGA.meth)]), ]
dim(data)
data[1:5,1:5]

dim(data)

# split data into numbers and probe info
data.only <- t(data[,-(1:4)])
data.info <- data[,1:4]
rm(data)


# setwd("/Users/12083/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/TCGA")
# getwd()


clinical.data <- read.delim("/mnt/ceph/fu_lab/TCGA/cbioportal/brca_tcga/data_bcr_clinical_data_patient.txt", skip = 4, header = TRUE)
clinical.data[1:5,1:5]
dim(clinical.data)

tmp <- clinical.data$`PATIENT_ID`[which(clinical.data$`ER_STATUS_BY_IHC` == "Positive")]
length(tmp)

table(clinical.data$`ER_STATUS_BY_IHC`)


for(i in 1:nrow(data.only)){

  print(i)

  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(rownames(data.only)[i]), '-'))
  split.ID

  #now use this ID to search in the clinical data
  new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3], sep = "")
  new.rowname

  #find the patient ID in the clinical data
  row.num <- which(as.character(new.rowname) == as.character(tmp))
  row.num

  if(length(row.num) > 0){

    age <- clinical.data$AGE[row.num]
    age

    race <- clinical.data$RACE[row.num]
    race

    final <- c(rownames(data.only)[i], age, race)
    write.table(t(final), file = "/mnt/ceph/kark6289/TCGA_analysis/names.pos.patient2.txt", sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE,quote=FALSE)
  }


}



neg <- fread("/mnt/ceph/kark6289/TCGA_analysis/names.neg.patient.txt", header =TRUE)

for(i in 1:nrow(neg)){
  print(i)

  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(neg[i]), '-'))
  split.ID

  #now merge the first 4 parts
  new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3],"-01", sep = "")
  new.rowname

  #remove this part and find age and race to use later
  new.rowname.2 <- paste(split.ID[1],".",split.ID[2],".",split.ID[3],".01", sep = "")
  new.rowname.2

  data <- cbind(neg[i], new.rowname,new.rowname.2)

  #save it to a new txt file
  write.table(data, file = "/mnt/ceph/kark6289/TCGA_analysis/split.names.neg.patient.txt", sep = "\t", row.names = FALSE,
              col.names = FALSE, append = TRUE, quote=FALSE)
}

pos <- fread("/mnt/ceph/kark6289/TCGA_analysis/names.pos.patient.txt", header =TRUE)

for(i in 1:nrow(pos)){

  print(i)

  #split the name in the TCGA data
  split.ID <- unlist(strsplit(as.character(pos[i]), '-'))
  split.ID

  #now merge the first 4 parts
  new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3],"-01", sep = "")
  new.rowname

  new.rowname.2 <- paste(split.ID[1],".",split.ID[2],".",split.ID[3],".01", sep = "")
  new.rowname.2

  data <- cbind(pos[i], new.rowname,new.rowname.2)

  #save it to a new txt file
  write.table(data, file = "/mnt/ceph/kark6289/TCGA_analysis/split.names.pos.patient.txt", sep = "\t", row.names = FALSE,
              col.names = FALSE, append = TRUE, quote=FALSE)
}






