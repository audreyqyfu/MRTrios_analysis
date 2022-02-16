
#Title: Simulation of M1 model data to compare between LOND and ADDIS

#Description: Constructing some simulations to compare between LOND and ADDIS 

#Date: 01-26-2021

#last updated: 02-02-2021


#=============================================================================================================


library(MRPC)

full_table <- function(M, N, signal){

#M <- 1000

#N <-  50

#signal <- 1

data_st <- list()

mrpc_lond_st <- list()

mrpc_addis_st <- list()

Truth_st <- list()


prerec_lond <- vector(length = M,
                      mode = 'list')

prerec_addis <- vector(length = M,
                       mode = 'list')

for (e in 1:M) {

  data <- SimulateData(N = N, 
                                       p = 0.45,
                                       'model4', 
                                       b0.1 = 0,
                                       b1.1 = signal, 
                                       b1.2 = signal,
                                       b1.3 = signal, 
                                       sd.1 = 1)

 # data <- data_1[sample(1:1000, N), ]
  
head(data)

# 'model1', 

  data_st[[e]] <- data
  
  
  n <- nrow (data)        # Number of rows
  V <- colnames(data)     # Column names
  
  # Calculate Pearson correlation
  suffStat_C <- list(C = cor(data),
                     n = n)
  
  # Infer the graph by MRPC
  mrpc_lond <- MRPC(data,
                    suffStat = suffStat_C,
                    GV = 1,
                    FDR = 0.05,
                    indepTest = 'gaussCItest',
                    labels = V,
                    FDRcontrol = 'LOND',
                    verbose = FALSE)
  
  mrpc_lond_st[[e]] <- mrpc_lond
  
  mrpc_addis <- MRPC(data,
                     suffStat = suffStat_C,
                     GV = 1,
                     FDR = 0.05,
                     indepTest = 'gaussCItest',
                     labels = V,
                     lambda = 0.5,
                     FDRcontrol = 'ADDIS',
                     verbose = FALSE)
  
  mrpc_addis_st[[e]] <- mrpc_addis
  
  #Adj_directed <- as(Truth, "matrix")
  
  library(help=MRPC)
  
  
  Truth <- MRPCtruth$M4   # Truth for model M4
  
  Truth_st[[e]] <- Truth
  
  #Adj_directed <- as(Truth, "matrix")
  
  prerec_lond[[e]] <- RecallPrecision(g1 = Truth,
                                      g2 = mrpc_lond@graph,
                                      GV = 1,
                                      includeGV = FALSE)
  

  
 # print(mrpc_lond@graph)
  prerec_addis[[e]] <- RecallPrecision(g1 =Truth,
                                       g2 = mrpc_addis@graph,
                                       GV = 1,
                                       includeGV = FALSE)
  
#print(mrpc_addis@graph)
}


  
  prerec <- matrix(c(mean(unlist(lapply(prerec_lond, `[`, 'Precision'))),
                         sd(unlist(lapply(prerec_lond, `[`, 'Precision'))),
                        mean(unlist(lapply(prerec_lond, `[`, 'Recall'))),
                         sd(unlist(lapply(prerec_lond, `[`, 'Recall'))),
                         mean(unlist(lapply(prerec_addis, `[`, 'Precision'))),
                         sd(unlist(lapply(prerec_addis, `[`, 'Precision'))),
                         mean(unlist(lapply(prerec_addis, `[`, 'Recall'))),
                         sd(unlist(lapply(prerec_addis, `[`, 'Recall')))),
                       byrow = TRUE,
                       nrow = 2)
  


colnames(prerec)<- colnames(prerec, do.NULL = FALSE)
colnames(prerec)<-c("Precision_mean","Precision_SD","Recall_mean", "Recall_SD") 
rownames(prerec)<- rownames(prerec, do.NULL = FALSE)
rownames(prerec)<- c("lond","Addis")

return(prerec) 

}

#generate the data for different M, N, signals 

#i for sample sizes 
#j for signals 

#M1
for (i in c(50, 200, 500, 1000 )) {
  
  for (j in c(0.2, 0.5, 1)) {
    
  print(i)
    
    print(j)
    
    print(full_table(1000,i, j)[1, 3:4]) 
  
}

}



full_table(10,50, 0.5) [1, 3:4]
full_table(1000,50, 1) [1, 3:4]



#save(prerec, file = '~/prerec_M1.RData')
 
#========================================================================================================
#Digging in the difference 

perc_lond_1 <- unlist(lapply(prerec_lond, `[`, 'Precision'))
perc_addis_1 <- unlist(lapply(prerec_addis, `[`, 'Precision'))

reca_lond_1 <- unlist(lapply(prerec_lond, `[`, 'Recall'))
reca_addis_1 <- unlist(lapply(prerec_addis, `[`, 'Recall'))

#Number of differences
length(which(perc_lond_1 !=perc_addis_1 ))

length(which(reca_lond_1 !=reca_addis_1 ))

#Looking into the differences  
which(reca_lond_1 !=reca_addis_1 )

x <- which(perc_lond_1 != perc_addis_1)[[1]]

#looking into first different data set 

head(data_st[[x]])

mrpc_lond_st[[x]]
mrpc_addis_st[[x]]


#plotting the difference
par(mfrow=c(1,3))

plot(Truth_st[[x]])
plot(mrpc_lond_st[[x]])
plot(mrpc_addis_st[[x]])


prerec_lond[[x]]

prerec_addis[[x]]



which(perc_lond_1 != perc_addis_1)

length(which(reca_lond_1 !=reca_addis_1 ))





for (x in 1:length(which(reca_lond_1 !=reca_addis_1 )) ) {
 
  head(data_st[[x]])
  
  mrpc_lond_st[[x]]
  mrpc_addis_st[[x]]
  
  plot(mrpc_lond_st[[x]])
  plot(mrpc_addis_st[[x]])
  
  
 print(prerec_lond[[x]])
  
  prerec_addis[[x]]
   
}







