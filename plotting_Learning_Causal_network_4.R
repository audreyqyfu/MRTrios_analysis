
# Title: Plotting the simulations 

# Description: Plotting the simulation comparision of different models 

# Created by: Mohamed Megheib

# Date: 01-01-2021

# Last updated: 02-23-2021


#=========================================================================================================

# Reading the data

par(mfrow=c(1,1))

pdf("LOND_ADDIS.pdf")

library(readxl)
M1_simulation_lond_Addis_Pre <- read_excel("C:/Users/megheib/OneDrive - University of Idaho/Lenovo_Yoga_08_05_2020/IMCI/Audrey/MRPC_Evan/LOND_ADDIS_Plots/M1_simulation_lond_Addis_Pre.xlsx")

#M1_simulation_lond_Addis <- read_excel("M2_simulation_lond_Addis.xlsx")

M1_simulation_lond_Addis_Pre <- read_excel("M1_simulation_lond_Addis_Pre.xlsx")
M1_simulation_lond_Addis_Rec <- read_excel("M1_simulation_lond_Addis_Rec.xlsx")

Complex_simulation_lond_Addis_Pre <- read_excel("Complex_simulation_lond_Addis_Pre.xlsx")
Complex_simulation_lond_Addis_Rec <- read_excel("Complex_simulation_lond_Addis_Rec.xlsx")


M0_simulation_lond_Addis_Pre <- read_excel("M0_simulation_lond_Addis_Pre.xlsx")
M0_simulation_lond_Addis_Rec <- read_excel("M0_simulation_lond_Addis_Rec.xlsx")




df_M1_Pre <- data.frame(M1_simulation_lond_Addis_Pre)
df_M1_Rec <- data.frame(M1_simulation_lond_Addis_Rec)





#require(ggplot2)
library(ggplot2)


colors <- c("LOND/ADDIS_0.2" = "blue", "LOND/ADDIS_0.5" = "red", "LOND/ADDIS_1" = "orange")

#M1

#Precision 
ggplot(df_M1_Pre, aes(x = c(50, 200, 500, 1000))) +
  geom_line(aes(y = df_M1_Pre$LOND_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M1_Pre$LOND_0.2_U, ymin = df_M1_Pre$LOND_0.2_L, color = "LOND/ADDIS_0.2"),linetype=1 )+
  
  geom_line(aes(y = df_M1_Pre$LOND_0.5_mean,  color = "LOND/ADDIS_0.5"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M1_Pre$LOND_0.5_U, ymin = df_M1_Pre$LOND_0.5_L, color = "LOND/ADDIS_0.5"),linetype=1 )+
  
  geom_line(aes(y = df_M1_Pre$LOND_1_mean, color = "LOND/ADDIS_1"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M1_Pre$LOND_1_U, ymin = df_M1_Pre$LOND_1_L, color = "LOND/ADDIS_1"),linetype=1 )+
 
  ##############ADDIS
  
  geom_line(aes(y = df_M1_Pre$ADDIS_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M1_Pre$ADDIS_0.2_U, ymin = df_M1_Pre$ADDIS_0.2_L, color = "LOND/ADDIS_0.2"),linetype=2 )+
  
  geom_line(aes(y = df_M1_Pre$ADDIS_0.5_mean,color = "LOND/ADDIS_0.5"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M1_Pre$ADDIS_0.5_U, ymin = df_M1_Pre$ADDIS_0.5_L,color = "LOND/ADDIS_0.5"),linetype=2)+
  
  geom_line(aes(y = df_M1_Pre$ADDIS_1_mean,color = "LOND/ADDIS_1"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M1_Pre$ADDIS_1_U, ymin = df_M1_Pre$ADDIS_1_L,color = "LOND/ADDIS_1"),linetype=3 )+
  scale_color_manual(values = colors)+
  xlab('sample_size') +
  ylab('Pre')+
  ggtitle("Precision_M1_Line is LOND and dots is ADDIS")

#Recall

ggplot(df_M1_Rec, aes(x = c(50, 200, 500, 1000))) +
  geom_line(aes(y = df_M1_Rec$LOND_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M1_Rec$LOND_0.2_U, ymin = df_M1_Rec$LOND_0.2_L, color = "LOND/ADDIS_0.2"),linetype=1 )+
  
  geom_line(aes(y = df_M1_Rec$LOND_0.5_mean,  color = "LOND/ADDIS_0.5"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M1_Rec$LOND_0.5_U, ymin = df_M1_Rec$LOND_0.5_L, color = "LOND/ADDIS_0.5"),linetype=1 )+
  
  geom_line(aes(y = df_M1_Rec$LOND_1_mean, color = "LOND/ADDIS_1"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M1_Rec$LOND_1_U, ymin = df_M1_Rec$LOND_1_L, color = "LOND/ADDIS_1"),linetype=1 )+
  
  ##############ADDIS
  
  geom_line(aes(y = df_M1_Rec$ADDIS_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M1_Rec$ADDIS_0.2_U, ymin = df_M1_Rec$ADDIS_0.2_L, color = "LOND/ADDIS_0.2"),linetype=2 )+
  
  geom_line(aes(y = df_M1_Rec$ADDIS_0.5_mean,color = "LOND/ADDIS_0.5"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M1_Rec$ADDIS_0.5_U, ymin = df_M1_Rec$ADDIS_0.5_L,color = "LOND/ADDIS_0.5"),linetype=2)+
  
  geom_line(aes(y = df_M1_Rec$ADDIS_1_mean,color = "LOND/ADDIS_1"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M1_Rec$ADDIS_1_U, ymin = df_M1_Rec$ADDIS_1_L,color = "LOND/ADDIS_1"),linetype=3 )+
  scale_color_manual(values = colors)+
  xlab('sample_size') +
  ylab('Rec')+
  ggtitle("Recall_M1_Line is LOND and dots is ADDIS")


######################################Complex

#Complex

df_Complex_Pre <- data.frame(Complex_simulation_lond_Addis_Pre)
df_Complex_Rec <- data.frame(Complex_simulation_lond_Addis_Rec)



#Precision 


ggplot(df_Complex_Pre, aes(x = c(50, 200, 500, 1000))) +
  geom_line(aes(y = df_Complex_Pre$LOND_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_Complex_Pre$LOND_0.2_U, ymin = df_Complex_Pre$LOND_0.2_L, color = "LOND/ADDIS_0.2"),linetype=1 )+
  
  geom_line(aes(y = df_Complex_Pre$LOND_0.5_mean,  color = "LOND/ADDIS_0.5"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_Complex_Pre$LOND_0.5_U, ymin = df_Complex_Pre$LOND_0.5_L, color = "LOND/ADDIS_0.5"),linetype=1 )+
  
  geom_line(aes(y = df_Complex_Pre$LOND_1_mean, color = "LOND/ADDIS_1"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_Complex_Pre$LOND_1_U, ymin = df_Complex_Pre$LOND_1_L, color = "LOND/ADDIS_1"),linetype=1 )+
  
  ##############ADDIS
  
  geom_line(aes(y = df_Complex_Pre$ADDIS_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_Complex_Pre$ADDIS_0.2_U, ymin = df_Complex_Pre$ADDIS_0.2_L, color = "LOND/ADDIS_0.2"),linetype=2 )+
  
  geom_line(aes(y = df_Complex_Pre$ADDIS_0.5_mean,color = "LOND/ADDIS_0.5"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_Complex_Pre$ADDIS_0.5_U, ymin = df_Complex_Pre$ADDIS_0.5_L,color = "LOND/ADDIS_0.5"),linetype=2)+
  
  geom_line(aes(y = df_Complex_Pre$ADDIS_1_mean,color = "LOND/ADDIS_1"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_Complex_Pre$ADDIS_1_U, ymin = df_Complex_Pre$ADDIS_1_L,color = "LOND/ADDIS_1"),linetype=3 )+
  scale_color_manual(values = colors)+
  xlab('sample_size') +
  ylab('Pre')+
  ggtitle("Precision_Complex_Line is LOND and dots is ADDIS")

#Recall

ggplot(df_Complex_Rec, aes(x = c(50, 200, 500, 1000))) +
  geom_line(aes(y = df_Complex_Rec$LOND_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_Complex_Rec$LOND_0.2_U, ymin = df_Complex_Rec$LOND_0.2_L, color = "LOND/ADDIS_0.2"),linetype=1 )+
  
  geom_line(aes(y = df_Complex_Rec$LOND_0.5_mean,  color = "LOND/ADDIS_0.5"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_Complex_Rec$LOND_0.5_U, ymin = df_Complex_Rec$LOND_0.5_L, color = "LOND/ADDIS_0.5"),linetype=1 )+
  
  geom_line(aes(y = df_Complex_Rec$LOND_1_mean, color = "LOND/ADDIS_1"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_Complex_Rec$LOND_1_U, ymin = df_Complex_Rec$LOND_1_L, color = "LOND/ADDIS_1"),linetype=1 )+
  
  ##############ADDIS
  
  geom_line(aes(y = df_Complex_Rec$ADDIS_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_Complex_Rec$ADDIS_0.2_U, ymin = df_Complex_Rec$ADDIS_0.2_L, color = "LOND/ADDIS_0.2"),linetype=2 )+
  
  geom_line(aes(y = df_Complex_Rec$ADDIS_0.5_mean,color = "LOND/ADDIS_0.5"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_Complex_Rec$ADDIS_0.5_U, ymin = df_Complex_Rec$ADDIS_0.5_L,color = "LOND/ADDIS_0.5"),linetype=2)+
  
  geom_line(aes(y = df_Complex_Rec$ADDIS_1_mean,color = "LOND/ADDIS_1"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_Complex_Rec$ADDIS_1_U, ymin = df_Complex_Rec$ADDIS_1_L,color = "LOND/ADDIS_1"),linetype=3 )+
  scale_color_manual(values = colors)+
  xlab('sample_size') +
  ylab('Rec')+
  ggtitle("Recall_Complex_Line is LOND and dots is ADDIS")





######################################M0

#M0

df_M0_Pre <- data.frame(M0_simulation_lond_Addis_Pre)
df_M0_Rec <- data.frame(M0_simulation_lond_Addis_Rec)



#Precision 


ggplot(df_M0_Pre, aes(x = c(50, 200, 500, 1000))) +
  geom_line(aes(y = df_M0_Pre$LOND_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M0_Pre$LOND_0.2_U, ymin = df_M0_Pre$LOND_0.2_L, color = "LOND/ADDIS_0.2"),linetype=1 )+
  
  geom_line(aes(y = df_M0_Pre$LOND_0.5_mean,  color = "LOND/ADDIS_0.5"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M0_Pre$LOND_0.5_U, ymin = df_M0_Pre$LOND_0.5_L, color = "LOND/ADDIS_0.5"),linetype=1 )+
  
  geom_line(aes(y = df_M0_Pre$LOND_1_mean, color = "LOND/ADDIS_1"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M0_Pre$LOND_1_U, ymin = df_M0_Pre$LOND_1_L, color = "LOND/ADDIS_1"),linetype=1 )+
  
  ##############ADDIS
  
  geom_line(aes(y = df_M0_Pre$ADDIS_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M0_Pre$ADDIS_0.2_U, ymin = df_M0_Pre$ADDIS_0.2_L, color = "LOND/ADDIS_0.2"),linetype=2 )+
  
  geom_line(aes(y = df_M0_Pre$ADDIS_0.5_mean,color = "LOND/ADDIS_0.5"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M0_Pre$ADDIS_0.5_U, ymin = df_M0_Pre$ADDIS_0.5_L,color = "LOND/ADDIS_0.5"),linetype=2)+
  
  geom_line(aes(y = df_M0_Pre$ADDIS_1_mean,color = "LOND/ADDIS_1"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M0_Pre$ADDIS_1_U, ymin = df_M0_Pre$ADDIS_1_L,color = "LOND/ADDIS_1"),linetype=3 )+
  scale_color_manual(values = colors)+
  xlab('sample_size') +
  ylab('Pre')+
  ggtitle("Precision_M0_Line is LOND and dots is ADDIS")

#Recall

ggplot(df_M0_Rec, aes(x = c(50, 200, 500, 1000))) +
  geom_line(aes(y = df_M0_Rec$LOND_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M0_Rec$LOND_0.2_U, ymin = df_M0_Rec$LOND_0.2_L, color = "LOND/ADDIS_0.2"),linetype=1 )+
  
  geom_line(aes(y = df_M0_Rec$LOND_0.5_mean,  color = "LOND/ADDIS_0.5"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M0_Rec$LOND_0.5_U, ymin = df_M0_Rec$LOND_0.5_L, color = "LOND/ADDIS_0.5"),linetype=1 )+
  
  geom_line(aes(y = df_M0_Rec$LOND_1_mean, color = "LOND/ADDIS_1"),size = 2, linetype=1) +
  geom_errorbar(aes(ymax = df_M0_Rec$LOND_1_U, ymin = df_M0_Rec$LOND_1_L, color = "LOND/ADDIS_1"),linetype=1 )+
  
  ##############ADDIS
  
  geom_line(aes(y = df_M0_Rec$ADDIS_0.2_mean, color = "LOND/ADDIS_0.2"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M0_Rec$ADDIS_0.2_U, ymin = df_M0_Rec$ADDIS_0.2_L, color = "LOND/ADDIS_0.2"),linetype=2 )+
  
  geom_line(aes(y = df_M0_Rec$ADDIS_0.5_mean,color = "LOND/ADDIS_0.5"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M0_Rec$ADDIS_0.5_U, ymin = df_M0_Rec$ADDIS_0.5_L,color = "LOND/ADDIS_0.5"),linetype=2)+
  
  geom_line(aes(y = df_M0_Rec$ADDIS_1_mean,color = "LOND/ADDIS_1"),size = 2, linetype=3) +
  geom_errorbar(aes(ymax = df_M0_Rec$ADDIS_1_U, ymin = df_M0_Rec$ADDIS_1_L,color = "LOND/ADDIS_1"),linetype=3 )+
  scale_color_manual(values = colors)+
  xlab('sample_size') +
  ylab('Rec')+
  ggtitle("Recall_M0_Line is LOND and dots is ADDIS")




dev.off()







