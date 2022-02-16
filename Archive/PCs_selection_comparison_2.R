#Title: Comparing different PCs selection approaches

# Description: Comparing different PCs selection approaches by ER test status in order to identify the best approach  

# Created by Mohamed Megheib

#Date: 11-01-2021

#Last updated: 12-10-2021

#============================================================================================================

# Independent setup 
Indep.setup.TCGA.8 <- read.csv("C:/Users/molot/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/TCGA results/11-11-2021/Indep.setup.TCGA.8.csv", header=FALSE)

# Joint setup 
Joint.setup.TCGA.6 <- read.csv("C:/Users/molot/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/TCGA results/11-11-2021/Joint.setup.TCGA.6.csv", header=FALSE)


dim(Joint.setup.TCGA.6)

Joint.setup.TCGA.6[1:5, 1:5]

dim(Indep.setup.TCGA.8)

Indep.setup.TCGA.8[1:4, 1:5]

Indep.setup.TCGA.8.1 <- Indep.setup.TCGA.8[1:dim(Joint.setup.TCGA.6)[1], ]

dim(Indep.setup.TCGA.8.1)

plot(Joint.setup.TCGA.6[, 4],Indep.setup.TCGA.8.1[, 4],type="p",  main="Comparison between joint (average=18.1) and
     independent (average=17.7) setup ", xlab="Joint", ylab="Indep" )

abline(0,1, col="red")

mean(Joint.setup.TCGA.6[, 4])
mean(Indep.setup.TCGA.8.1[, 4])
#============================================Indep_1st versus Joint==========================================

# Reading the data 
Joint_Indep_Comparison_2 <- read.delim("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Joint_Indep_PCA/Joint_Indep_Comparison_3.txt")
dim(Joint_Indep_Comparison_2)
tail(Joint_Indep_Comparison_2)
names(Joint_Indep_Comparison_2)

# Total 
plot(Joint_Indep_Comparison_2$indep.total.3,Joint_Indep_Comparison_2$joint.total.3,type="p",  main="Comparison between total joint and
     1st indep setup ",xlab="Indep", ylab="Joint" )
abline(0,1, col="red")
mean(Joint_Indep_Comparison_2$indep.total.3)
mean(Joint_Indep_Comparison_2$joint.total.3)

# Positive 
plot(Joint_Indep_Comparison_2$indep.positive.3 ,Joint_Indep_Comparison_2$joint.positive.3,type="p",  main="Comparison between positive joint and
     1st indep setup ",xlab="Indep", ylab="Joint" )
abline(0,1, col="red")
mean(Joint_Indep_Comparison_2$indep.positive.3)
mean(Joint_Indep_Comparison_2$joint.positive.3)

#Negative 
plot(Joint_Indep_Comparison_2$indep.negative.3,Joint_Indep_Comparison_2$joint.negative.3,type="p",  main="Comparison between negative joint and
     1st indep setup ",xlab="Indep", ylab="Joint" )
abline(0,1, col="red")
mean(Joint_Indep_Comparison_2$indep.negative.3)
mean(Joint_Indep_Comparison_2$joint.negative.3)

length(Joint_Indep_Comparison_2$indep.negative.3)
mean(Joint_Indep_Comparison_2$joint.negative.3)

#================================================Indep_2nd versus Joint============================================

#Reading the data 
Indep.setup.TCGA.8.2nd <- read.delim("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Joint_Indep_PCA/Indep_2nd/Indep.setup.TCGA.8.2nd.txt", header=FALSE)
dim(Indep.setup.TCGA.8.2nd)
Indep.setup.TCGA.8.2nd[1:3, 1:4]

#Total
plot(Indep.setup.TCGA.8.2nd[1:dim(Joint_Indep_Comparison_2)[1], ]$V4,Joint_Indep_Comparison_2$joint.total.3,type="p",  main="Comparison between total joint and
     2nd indep setup ",xlab="Indep", ylab="Joint" )
abline(0,1, col="red")
mean(Indep.setup.TCGA.8.2nd[1:dim(Joint_Indep_Comparison_2)[1], ]$V4)
mean(Joint_Indep_Comparison_2$joint.total.3)

#Positive 
Indep.setup.TCGA.8.2nd.Positive <- read.delim("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Joint_Indep_PCA/Indep_2nd/Indep.setup.TCGA.8.2nd.Positive.txt", header=FALSE)
dim(Indep.setup.TCGA.8.2nd.Positive)

plot(Indep.setup.TCGA.8.2nd.Positive[1:dim(Joint_Indep_Comparison_2)[1], ]$V4 ,Joint_Indep_Comparison_2$joint.positive.3,type="p",  main="Comparison between positive joint and
     2st indep setup ",xlab="Indep", ylab="Joint" )
abline(0,1, col="red")
mean(Indep.setup.TCGA.8.2nd.Positive[1:dim(Joint_Indep_Comparison_2)[1], ]$V4)
mean(Joint_Indep_Comparison_2$joint.positive.3)


#Negative
Indep.setup.TCGA.8.2nd.Negative <- read.delim("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/Joint_Indep_PCA/Indep_2nd/Indep.setup.TCGA.8.2nd.Negative.txt", header=FALSE)
dim(Indep.setup.TCGA.8.2nd.Negative)
plot(Indep.setup.TCGA.8.2nd.Negative[1:dim(Joint_Indep_Comparison_2)[1], ]$V4,Joint_Indep_Comparison_2$joint.negative.3,type="p",  main="Comparison between negative joint and
     2st indep setup ",xlab="Indep", ylab="Joint" )
abline(0,1, col="red")
mean(Indep.setup.TCGA.8.2nd.Negative[1:dim(Joint_Indep_Comparison_2)[1], ]$V4)
mean(Joint_Indep_Comparison_2$joint.negative.3)

#===================================================================Indep_1st versus Indep_2nd===============================


plot(Indep.setup.TCGA.8.2nd[1:dim(Joint_Indep_Comparison_2)[1], ]$V4,Joint_Indep_Comparison_2$indep.total.3,type="p",  main="Comparison between total 1st indep setup and
     2nd indep setup ",xlab="Indep_2nd", ylab="Indep_1st" )
abline(0,1, col="red")



plot(Indep.setup.TCGA.8.2nd.Positive[1:dim(Joint_Indep_Comparison_2)[1], ]$V4,Joint_Indep_Comparison_2$joint.positive.3,type="p",  main="Comparison between positive 1st indep setup and
     2nd indep setup ",xlab="Indep_2nd", ylab="Indep_1st" )
abline(0,1, col="red")


plot(Indep.setup.TCGA.8.2nd.Negative[1:dim(Joint_Indep_Comparison_2)[1], ]$V4,Joint_Indep_Comparison_2$joint.negative.3,type="p",  main="Comparison between total 1st indep setup and
     2nd indep setup ",xlab="Indep_2nd", ylab="Indep_1st" )
abline(0,1, col="red")





