#PROJECT: Investigating Causal Network Inference for Transcription and Methylation Relationship in Breast Cancer

###09.04.2024###
#set working directory
setwd("C:/Users/02kri/OneDrive - Wayne State University")

##installed data table package using this code:
#install.packages('data.table')

##Activating trios data with data.frame
library('data.table')

#Loading in causal model inference of ER+(code reference: methylation_data)
methylation_data <-read.delim('model.trio.MRGN.all.posER.reclassify2.txt')

#Loading in list of TRIOS consisting of gene expression (GE), methylation and   (code reference: trios)
trios <- data.frame(fread('trio.final.protein.coding.txt'))

##Data Exercises to manipulate data (SAVE just in case)
#head(trios)
#trios[2,]
#trios[2,3]
#head(methylation_data)
#trios[270001,]
#trios[270002,]
#trios[270001:270010,]
#methylation_data$trio_row[1:10]
#unique_genes <- trios[methylation_data$trio_row[1:10],]
#unique_genes$Gene.name
#unique(unique_genes$Gene.name)
#length(unique(unique_genes$Gene.name))



###09.07.2024###
#Finding overall unique genes within ER+ Data
positive_genes <- trios[methylation_data$trio_row,]
results <-unique(positive_genes$Gene.name)
length(results) # = 16331
length(unique(positive_genes$gene.row)) #16328? (off by 3)
##OBSERVATION: gene name & gene row columns in ER+ data do not align in the # of unique genes (possibly merging genes together?)

#Loading in causal model inference of ER- (code reference: second_methylation_data)
second_methylation_data <- read.delim('model.trio.MRGN.all.negER.reclassify2.txt')

#Finding overall unique genes within ER- Data
negative_genes <-trios[second_methylation_data$trio_row,]
other_results <- unique(negative_genes$Gene.name) 
length(other_results) # = 16282
length(unique(negative_genes$gene.row)) #16279? (also off by 3)
##OBSERVATION: same pattern goes for ER- data



###09.17.2024###
#Loading in raw data for TRIOS (use 'fread' to load in large & heavily filtered datasets)
CNA_rawdata <- fread('data_CNA.txt') #CNA
GE_rawdata <- fread('data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt') #geneexpression
meth_rawdata <- fread('split.names.TCGA.meth.logit.txt') #methylation

#Loading in all intermediate files
##ER+ intermediate files (PCA files will have genes as rows instead of columns)
ERpos_patients <- read.delim('names.pos.patient2.txt', sep = ',')#ER+ patient list
pos_GE_PCA <- read.delim('PCA.gene.exp.posER.txt', row.names = 1) #gene expression ER+ PCA scores
pos_meth_PCA <- read.delim('PCA.meth.posER.txt', row.names = 1) #methylation ER+ PCA scores

##ER- intermediate files
#ERneg_patients <- read.delim('names.neg.patient2.txt', sep = ',') #ER- patient list
#neg_GE_PCA <- read.delim('PCA.gene.exp.negER.txt', row.names = 1) #gene expression ER- PCA scores
#neg_meth_PCA <-read.delim('PCA.meth.negER.txt', row.names = 1) #methylation ER- PCA scores


###09.18.2024###
#GOAL: Generate pairwise scatter plots to visualize causal models (code reference: "methylation_data") 
#PURPOSE: Must investigate extracted Trios with genomic data at individual levels

#Avoiding errors of "list" by turning into dataframe
GE_rawdata_tmp <- as.data.frame(GE_rawdata) #made with Dr. Fu
ME_rawdata_tmp <- as.data.frame(meth_rawdata)
CNA_rawdata_tmp <- as.data.frame(CNA_rawdata)

#Extracting the row names of ER+ methylation PCA intermediate file 
#this is the patients' IDs
first_column <- rownames(pos_meth_PCA)

#Matching columns of GE_raw data to extracted row names of ER+ methylation PCA intermediate file
final_matching <- match(first_column, colnames(GE_rawdata))
#final_matching
GE_ERpos <- GE_rawdata_tmp[,final_matching]



###09.23.2024###
#Extracting the row names of ER+ GE PCA intermediate file
#this is the patients' IDs
column <- rownames(pos_GE_PCA)

#Matching columns of GE_raw data to extracted row names of ER+ GE PCA intermediate file 
matching <- match(column, colnames(GE_rawdata))
#matching
ER_pos <- GE_rawdata_tmp[,matching]

##OBSERVATION: using the GE and methylation PCA intermediate files will give you the same result
#Code below proves this point:
first_column == column  #all elements are TRUE
##CONCLUSION: do not have to match all raw datas to both GE and methylation PCA files



###09.25.2024###
#NOTE: Sticking with "first column"/pos_meth_PCA when extracting row names (code: first_column <- rownames(pos_meth_PCA))
#Matching columns of meth_rawdata to extracted row names of ER+ methylation PCA intermediate file
more_matching <- match(first_column, colnames(meth_rawdata))
more_ER_pos <- ME_rawdata_tmp[, more_matching] 
  

###10.02.2024###
#RECAP: CNA data was fincky with matching column to PCA intermediate file
##Previous code: 
moree_matchingg <-match(first_column, colnames(CNA_rawdata))

##Why it did not work?:
#Check "moree_matchingg" vector with code: moree_matchingg[1:10]
#Conclusion: matching function is not a fan of "NA" and dimensions of GE_rawdata & CNA rawdata are different

##Question: How many NAs/Which NAs do we remove from CNA_rawdata?
#Code: sum(is.na(moree_matchingg))  = 8

#Revising first_column for CNA
revised_patients <- first_column[ -which (is.na (moree_matchingg))] #basically keeps "first_column" - all the NAs

#Matching extracted row names of ER- methylation PCA intermediate file to columns of CNA_rawdata
moree_matchingg <- match(revised_patients, colnames(CNA_rawdata))
positive_CNA <- CNA_rawdata_tmp[,moree_matchingg]


#Now the number of patients is different between each raw data set (572 variables v. 564 variables). 
#We must eliminate the same eight patients for all matching scenarios
###10.09.2024###
#NEW matching with GE data using "revised_patients"
matching <-match(revised_patients, colnames(GE_rawdata))
positive_GE<- GE_rawdata_tmp[,matching]

#NEW matching with meth raw data using "revised_patients"
more_matching <- match(revised_patients, colnames(meth_rawdata))
positive_meth <- ME_rawdata_tmp[, more_matching] 

##Using methylation_data: looking at all trios between mediation models

#For m1.1(C -> E -> M):
#example: the 12th row of methylation_data contains this inferred model
#code: methylation_data[12,]
#running the code, look at the "trio_row" column number: 270016
#now use the code: trios[270016,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
plot(unlist(positive_GE[9100,]), unlist (positive_meth[442663,]), main = 'm1.1 Gene Expression v. Methylation', xlab= 'Gene Expression', ylab = 'Methylation', col="blue", pch=16)
my.cor.m11 <- cor(unlist(positive_GE[9100,]), unlist (positive_meth[442663,]), use = 'pairwise.complete.obs')
legend ('topright', legend = paste('r = ', round(my.cor.m11, digits = 3)))

plot(unlist (positive_CNA[15092,]), unlist(positive_GE[9100,]), main ='m1.1 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alternation', ylab = 'Gene Expression', col = 'blue', pch = 16)
my.cor.m110<- cor(unlist (positive_CNA[15092,]), unlist(positive_GE[9100,]))
legend ('topright', legend = paste('r = ', round(my.cor.m110, digits = 3)))

plot(unlist (positive_CNA[15092,]), unlist(positive_meth[442663,]), main = 'm1.1 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alternation', ylab = 'Methylation', col = 'blue', pch = 16)
my.cor.m1100<- cor(unlist (positive_CNA[15092,]), unlist(positive_meth[442663,]), use = 'pairwise.complete.obs')
legend ('topright', legend = paste('r = ', round(my.cor.m1100, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = -0.561 
#CNA v. GE -> r = 0.217 
#CNA v. Methylation -> r = -0.132



###10.12.2024###
#For m1.2 (C -> M -> E):
#example: the 19th row of methylation_data contains this inferred model
#code: methylation_data[19,]
#now use the code: trios[270023,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
plot(unlist(positive_GE[5411,]), unlist (positive_meth[222723,]), main = 'Initial m1.2 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'deeppink', pch = 16)
my.cor.m12<- cor(unlist(positive_GE[5411,]), unlist (positive_meth[222723,]))
legend ('topright', legend = paste('r = ', round(my.cor.m12, digits = 3)))

plot(unlist (positive_CNA[20630,]), unlist(positive_GE[5411,]), main = 'Initial m1.2 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'deeppink', pch = 16)
my.cor.m120 <- cor(unlist (positive_CNA[20630,]), unlist(positive_GE[5411,]))
legend ('topright', legend = paste('r = ', round(my.cor.m120, digits = 3)))

plot(unlist (positive_CNA[20630,]), unlist(positive_meth[222723,]), main = 'Initial m1.2 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'deeppink', pch = 16)
my.cor.m1200 <- cor(unlist (positive_CNA[20630,]), unlist(positive_meth[222723,]))
legend ('topright', legend = paste('r = ', round(my.cor.m1200, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = 0.317 
#CNA v. GE -> r = -0.116 
#CNA v. Methylation -> r = 0.018 (essentially no correlation)


#For m0.1 (C → E; no relationship  between E and M):
#example: the 1st row of methylation_data contains this inferred model
#code: methylation_data[1,]
#now use the code: trios[270001,]
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
plot(unlist(positive_GE[1197,]), unlist (positive_meth[462699,]), main = 'm0.1 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'purple', pch = 16)
my.cor.m01 <- cor(unlist(positive_GE[1197,]), unlist (positive_meth[462699,]))
legend ('topright', legend = paste('r = ', round(my.cor.m01, digits = 3)))

plot(unlist(positive_CNA[264,]), unlist(positive_GE[1197,]), main = 'm0.1 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'purple', pch = 16)
my.cor.m010 <- cor(unlist(positive_CNA[264,]), unlist(positive_GE[1197,]))
legend ('topright', legend = paste('r = ', round(my.cor.m010, digits = 3)))

plot(unlist(positive_CNA[264,]), unlist(positive_meth[462699,]), main = 'm0.1 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'purple', pch = 16)
my.cor.m0100 <- cor(unlist (positive_CNA[264,]), unlist(positive_meth[462699,]))
legend ('topright', legend = paste('r = ', round(my.cor.m0100, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = 0.118 
#CNA v. GE -> r = 0.327 
#CNA v. Methylation -> r = 0.120 


#For m0.2 (C → M; no relationship  between E and M):
#example: the 31st row of methylation_data contains this inferred model
#code: methylation_data[31,]
#now use the code: trios[270035,] (NOTE: uses same GE and CNA value as m1.2)
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
plot(unlist(positive_GE[5411,]), unlist (positive_meth[424598,]), main = 'Initial m0.2 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'salmon', pch = 16)
my.cor.m02 <- cor(unlist(positive_GE[5411,]), unlist (positive_meth[424598,]))
legend ('topright', legend = paste('r = ', round(my.cor.m02, digits = 3)))

plot(unlist (positive_CNA[20630,]), unlist(positive_GE[5411,]), main = 'Initial m0.2 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'salmon', pch = 16)
my.cor.m020 <- cor(unlist (positive_CNA[20630,]), unlist(positive_GE[5411,]))
legend ('topright', legend = paste('r = ', round(my.cor.m020, digits = 3)))

plot(unlist(positive_CNA[20630,]), unlist(positive_meth[424598,]), main = 'Initial m0.2 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'salmon', pch = 16)
my.cor.m0200 <- cor(unlist (positive_CNA[20630,]), unlist(positive_meth[424598,]))
legend ('topright', legend = paste('r = ', round(my.cor.m0200, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed 
#GE v. Methylation -> r = 0.084 (essentially no correlation)
#CNA v. GE -> r = -0.116 
#CNA v. Methylation -> r = -0.263 


#For M3 (E <- C -> M):
#example: the 2nd row of methylation_data contains this inferred model
#code: methylation_data[2,]
#now use the code: trios[270002,] (NOTE: uses same GE value as m1.1)
#meth.row, cna.row and gene.row should contain numbers for you to use in scatterplots!
#make sure these values are not in "lists"
plot(unlist(positive_GE[9100,]), unlist (positive_meth[59380,]), main = 'Initial M3 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'red', pch = 16)
my.cor.m3 <- cor(unlist(positive_GE[9100,]), unlist (positive_meth[59380,]))
legend ('topright', legend = paste('r = ', round(my.cor.m3, digits = 3)))

plot(unlist(positive_CNA[15092,]), unlist(positive_GE[9100,]), main = 'Initial M3 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'red', pch = 16)
my.cor.m30 <- cor(unlist(positive_CNA[15092,]), unlist(positive_GE[9100,]))
legend ('topright', legend = paste('r = ', round(my.cor.m30, digits = 3)))

plot(unlist (positive_CNA[15092,]), unlist(positive_meth[59380,]), main = 'Initial M3 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'red', pch = 16)
my.cor.m300 <- cor(unlist (positive_CNA[15092,]), unlist(positive_meth[59380,]))
legend ('topright', legend = paste('r = ', round(my.cor.m300, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = 0.011 (essentially no correlation)
#CNA v. GE -> r = 0.217 
#CNA v. Methylation -> r = 0.092 (essentially no correlation)


###10.16.2024/10.20.2024###
##New m1.2 since EFNA2 is an outlier:
#New method of finding inferred models
#CODE: which(methylation_data$Inferred.Model3 == 'M1.2' & methylation_data$Inferred.Model2 == 'M1.2'& methylation_data$Inferred.Model == 'M1.2')
#code: methylation_data[632,]
#now use the code:trios[270643,]
plot(unlist(positive_GE[194,]), unlist(positive_meth[195911,]), main = 'm1.2 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col='deeppink', pch=16)
new.cor.m12 <- cor(unlist(positive_GE[194,]), unlist(positive_meth[195911,]))
legend ('topright', legend = paste('r = ', round(new.cor.m12, digits = 3)))

plot(unlist(positive_CNA[23833,]), unlist(positive_GE[194,]), main = 'm1.2 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col='deeppink', pch=16)
new.cor.m120 <-cor(unlist(positive_CNA[23833,]), unlist(positive_GE[194,]))
legend ('topright', legend = paste('r = ', round(new.cor.m120, digits = 3)))

plot(unlist(positive_CNA[23833,]), unlist(positive_meth[195911,]), main = 'm1.2 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col='deeppink', pch=16)
new.cor.m1200 <- cor(unlist(positive_CNA[23833,]), unlist(positive_meth[195911,]))
legend ('topright', legend = paste('r = ', round(new.cor.m1200, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = -0.205 
#CNA v. GE -> r =  0.144 
#CNA v. Methylation -> r = 0.088 (essentially no correlation)
#?? Expectation is that CNA v. GE has a stronger correlation (perhaps involves co-founders?)


#For new m0.2 since EFNA2 and others such as lines 249 & 396 is an outlier (SAME as m1.2)
#CODE: which(methylation_data$Inferred.Model3 == 'M0.2' & methylation_data$Inferred.Model2 == 'M0.2'& methylation_data$Inferred.Model == 'M0.2')
#code: methylation_data[1027,]
#now use the code:trios[271038,]
plot(unlist(positive_GE[952,]), unlist (positive_meth[261661,]), main = 'm0.2 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'salmon', pch = 16)
new.cor.m02 <- cor(unlist(positive_GE[952,]), unlist (positive_meth[261661,]))
legend ('topright', legend = paste('r = ', round(new.cor.m02, digits = 3)))

plot(unlist(positive_CNA[23750,]), unlist (positive_GE[952,]), main = 'm0.2 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'salmon', pch = 16)
new.cor.m020 <- cor(unlist(positive_CNA[23750,]), unlist (positive_GE[952,]))
legend ('topright', legend = paste('r = ', round(new.cor.m020, digits = 3)))

plot(unlist(positive_CNA[23750,]), unlist (positive_meth[261661,]), main = 'm0.2 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'salmon', pch = 16)
new.cor.m0200 <- cor(unlist(positive_CNA[23750,]), unlist (positive_meth[261661,]))
legend ('topright', legend = paste('r = ', round(new.cor.m0200, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed 
#GE v. Methylation -> r =  0.003 (essentially no correlation)
#CNA v. GE -> r = 0.018 (essentially no correlation)
#CNA v. Methylation -> r =  0.109 


#For m2.1 (C -> E <- M):
#CODE: which(methylation_data$Inferred.Model3 == 'M2.1' & methylation_data$Inferred.Model2 == 'M2.1'& methylation_data$Inferred.Model == 'M2.1')
#code: methylation_data[499,]
#now use the code:trios[270510,]
plot(unlist(positive_GE[3194,]), unlist (positive_meth[481639,]), main = 'm2.1 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'seagreen', pch = 16)
my.cor.m21 <- cor(unlist(positive_GE[3194,]), unlist (positive_meth[481639,]))
legend ('topright', legend = paste('r = ', round(my.cor.m21, digits = 3)))

plot(unlist(positive_CNA[13439,]), unlist (positive_GE[3194,]), main = 'm2.1 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'seagreen', pch = 16)
my.cor.m210 <- cor(unlist(positive_CNA[13439,]), unlist (positive_GE[3194,]))
legend ('topright', legend = paste('r = ', round(my.cor.m210, digits = 3)))

plot(unlist(positive_CNA[13439,]), unlist (positive_meth[481639,]), main = 'm2.1 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'seagreen', pch = 16)
my.cor.m2100 <-cor(unlist(positive_CNA[13439,]), unlist (positive_meth[481639,]))
legend ('topright', legend = paste('r = ', round(my.cor.m2100, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = -0.010 (essentially no correlation)
#CNA v. GE -> r = 0.379 
#CNA v. Methylation -> r = -0.042 (essentially no correlation)


#For m2.2 (C -> M <- E):
#CODE: which(methylation_data$Inferred.Model3 == 'M2.2' & methylation_data$Inferred.Model2 == 'M2.2'& methylation_data$Inferred.Model == 'M2.2')
#code: methylation_data[1116,]
#now use the code:trios[271127,]
plot(unlist(positive_GE[18832,]), unlist (positive_meth[422513,]), main = 'm2.2 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'green', pch = 16)
#cor(unlist(positive_GE[18832,]), unlist (positive_meth[422513,]))
legend('topright', legend = 'r = 0.043')

plot(unlist(positive_CNA[16315,]), unlist (positive_GE[18832,]), main = 'm2.2 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'green', pch = 16)
#cor(unlist(positive_CNA[16315,]), unlist (positive_GE[18832,]))
legend('topright', legend = 'r = 0.078')

plot(unlist(positive_CNA[16315,]), unlist (positive_meth[422513,]), main = 'm2.2 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'green', pch = 16)
#cor(unlist(positive_CNA[16315,]), unlist (positive_meth[422513,]))
legend('topright', legend = 'r = 0.221')

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = 0.043 (essentially no correlation)
#CNA v. GE -> r = 0.078 (essentially no correlation)
#CNA v. Methylation -> r = 0.221 


#For M4 (C -> M <- E):
#CODE: which(methylation_data$Inferred.Model3 == 'M2.2' & methylation_data$Inferred.Model2 == 'M2.2'& methylation_data$Inferred.Model == 'M2.2')
#code: methylation_data[785,]
#now use the code:trios[270796,]
plot(unlist(positive_GE[9372,]), unlist (positive_meth[483273,]), main = 'M4 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'orange', pch = 16)
my.cor.m4 <- cor(unlist(positive_GE[9372,]), unlist (positive_meth[483273,]), use = 'pairwise.complete.obs')
legend ('topright', legend = paste('r = ', round(my.cor.m4, digits = 3)))

plot(unlist(positive_CNA[16444,]), unlist (positive_GE[9372,]), main = 'M4 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'orange', pch = 16)
my.cor.m40 <- cor(unlist(positive_CNA[16444,]), unlist (positive_GE[9372,]))
legend ('topright', legend = paste('r = ', round(my.cor.m40, digits = 3)))

plot(unlist(positive_CNA[16444,]), unlist (positive_meth[483273,]), main = 'M4 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'orange', pch = 16)
my.cor.m400 <- cor(unlist(positive_CNA[16444,]), unlist (positive_meth[483273,]), use = 'pairwise.complete.obs')
legend ('topright', legend = paste('r = ', round(my.cor.m400, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = - 0.107 
#CNA v. GE -> r = 0.398 
#CNA v. Methylation -> r = -0.141 



###10.28.2024###
#For new M3:
#CODE: which(methylation_data$Inferred.Model3 == 'M3' & methylation_data$Inferred.Model2 == 'M3'& methylation_data$Inferred.Model == 'M3')
#code: methylation_data[804,]
#now use the code: trios[270815,]
plot(unlist(positive_GE[12246,]), unlist (positive_meth[356560,]), main = 'M3 Gene Expression v. Methylation', xlab = 'Gene Expression', ylab = 'Methylation', col = 'red', pch = 16)
new.cor.m3 <- cor(unlist(positive_GE[12246,]), unlist (positive_meth[356560,]))
legend ('topright', legend = paste('r = ', round(new.cor.m3, digits = 3)))

plot(unlist(positive_CNA[21113,]), unlist(positive_GE[12246,]), main = 'M3 Copy Number Alteration v. Gene Expression', xlab = 'Copy Number Alteration', ylab = 'Gene Expression', col = 'red', pch = 16)
new.cor.m30 <- cor(unlist(positive_CNA[21113,]), unlist(positive_GE[12246,]))
legend ('topright', legend = paste('r = ', round(new.cor.m30, digits = 3)))

plot(unlist(positive_CNA[21113,]), unlist(positive_meth[356560,]), main = 'M3 Copy Number Alteration v. Methylation', xlab = 'Copy Number Alteration', ylab = 'Methylation', col = 'red', pch = 16)
new.cor.m300 <-cor(unlist(positive_CNA[21113,]), unlist(positive_meth[356560,]))
legend ('topright', legend = paste('r = ', round(new.cor.m300, digits = 3)))

##Correlation (look at correlation as an absolute value): 
##+/- only shows what way the correlation is graphed
#GE v. Methylation -> r = -0.227 
#CNA v. GE -> r = 0.335 
#CNA v. Methylation -> r = -0.233 



###NOT ABLE TO APPLY TO ER- PATIENTS; REVISE AND APPLY THE SAME METHODS ABOVE...###

##Extracting the row names of ER- methylation PCA intermediate file
#second_column <- rownames(neg_meth_PCA)

##Matching extracted row names of ER- methylation PCA intermediate file to columns of GE_rawdata
#NOTE: Sticking with "second column"/neg_meth_PCA when extracting row names (code: second_column <- rownames(neg_meth_PCA))
#other_matching <- match(second_column, colnames(GE_rawdata))
#GE_ERneg <- GE_rawdata_tmp[ ,other_matching] 

##Matching extracted row names of ER- methylation PCA intermediate file to columns of meth_rawdata
#second_matching <- match(second_column, colnames(meth_rawdata))
#more_ER_neg <- ME_rawdata_tmp[, second_matching] 

##Matching extracted row names of ER- methylation PCA intermediate file to columns of CNA_rawdata
#secondd_matchingg <- match(second_column, colnames(CNA_rawdata))
#moree_neg <- CNA_rawdata_tmp[ , secondd_matchingg] #CODE ALSO NOT WORKING


