
########## histograms ###########

pdf(file = "/mnt/ceph/kark6289/TCGA_analysis/newhist/unscaled/TCGA.pdf", onefile=TRUE)

for(i in 1:100) {
  
  print(i)
  
  data.nona <- log(data[,5:899]/(1-data[,5:899]))
  
  hist(as.numeric(data.nona), xlim = c(-2,3.5), main = paste("Histogram of (TCGA) ",data$Row.names[i], sep = ""), xlab = "Methylation")
  
}

dev.off()

#################################################

setwd("/mnt/ceph/megheib/Cancer_Genomics/TCGA_data_results")
getwd()

load("TCGA.meth.RData")
dim(TCGA.meth)

data <- TCGA.meth[rowSums(is.na(TCGA.meth[,5:899])) != ncol(TCGA.meth[,5:899]), ]

dim(data)

# split data into numbers and probe info
data.only <- t(data[,-(1:4)])
data.info <- data[,1:4]
rm(data)


setwd("/Users/12083/OneDrive - University of Idaho/DOCUMENTS/AUDREY FU LAB/TCGA")
# getwd()


library(data.table)
clinical.data <- fread("brca_tcga_clinical_data.txt")
clinical.data[1:5,1:5]
dim(clinical.data)

library(tidyr)
tmp <- clinical.data %>% drop_na(`ER Status By IHC`)
tmp[1:5,1:5]
dim(tmp)

# rows <- grep("TCGA-A7-A26I",rownames(data.only))
# rows[1:5]

for(i in 1:nrow(data.only)){
  
#split the name in the TCGA data  
split.ID <- unlist(strsplit(as.character(rownames(data.only)[i]), '-'))
split.ID

#now use this ID to search in the clinical data
new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3], sep = "")
new.rowname

#find the patient ID in the clinical data
row.num <- which(as.character(new.rowname) == as.character(tmp$`Patient ID`))
row.num

if(length(row.num)>0){
  
if(tmp$`ER Status By IHC`[row.num] == "Positive"){
  
  ER.pos <- c(new.rowname,data.only[i,],tmp$`ER Status By IHC`[row.num])
  
  write.table(t(ER.pos), file = "/mnt/ceph/kark6289/TCGA_analysis/ERpos.csv", sep = ",", row.names = FALSE,
              col.names = FALSE, append = TRUE,quote=FALSE)
}else{
  
  ER.neg <- c(new.rowname,data.only[i,],tmp$`ER Status By IHC`[row.num])
  write.table(t(ER.neg), file = "/mnt/ceph/kark6289/TCGA_analysis/ERneg.csv", sep = ",", row.names = FALSE,
              col.names = FALSE, append = TRUE,quote=FALSE)

}
}
}

#######################################

setwd("/mnt/ceph/megheib/Cancer_Genomics/TCGA_data_results")
getwd()

load("TCGA.meth.RData")
dim(TCGA.meth)

data <- TCGA.meth[rowSums(is.na(TCGA.meth[,5:899])) != ncol(TCGA.meth[,5:899]), ]

dim(data)
data[1:5,1:5]

# split data into numbers and probe info
data.only <- t(data[,-(1:4)])
data.info <- data[,1:4]
rm(data)

final.col.num = NULL

#split individual ID and match
for(i in 1:nrow(data.only)){
  
split.ID <- unlist(strsplit(as.character(rownames(data.only)[i]), '-'))
split.ID

#now use this ID to search in the clinical data
new.rowname <- paste(split.ID[1],"-",split.ID[2],"-",split.ID[3],"-","01", sep = "")
new.rowname

#find the patient ID in the clinical data
col.num <- which(as.character(new.rowname) == as.character(colnames(cna)))
col.num

final.col.num <- c(final.col.num,col.num)

}

#match gene names
tmp <- na.omit(match(data$Gene_Symbol,cna$Hugo_Symbol))


# res <- cna[tmp,]
# final.res <- res[,final.col.num]

write.table(final, file = "/mnt/ceph/kark6289/TCGA_analysis/split.names.TCGA.meth.logit.txt", sep = ",", row.names = FALSE,
            col.names = TRUE, append = TRUE,quote=FALSE)




