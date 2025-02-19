
# Read the datasets from saved ones
M1.2_neg <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M1.2_neg_extract.txt",sep = "\t")
M1.2_pos <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M1.2_pos_extract.txt",sep = "\t")
M1.1_neg <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M1.1_neg_extract.txt",sep = "\t")
M1.1_pos <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M1.1_pos_extract.txt",sep = "\t")

M0.2_neg <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M0.2_neg_extract.txt",sep = "\t")
M0.2_pos <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M0.2_pos_extract.txt",sep = "\t")
M0.1_neg <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M0.1_neg_extract.txt",sep = "\t")
M0.1_pos <- read.delim("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Data_Result_Methyl_LZ/M0.1_pos_extract.txt",sep = "\t")


##########################################################################################################################
######  Distribution of the log10 (distance) between the nearby CpG island and the CpG in the methylation probe ##########
##########################################################################################################################

setwd("/Users/lianzuo/LZ/WayneSU/Methyl_Analysis_LZ/Plot_Result_Methyl_LZ")
library(ggplot2)
library(patchwork)
library(cowplot)

#hist_M0.1_pos_log_distance=ggplot(M0.1_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() , aes(x = log10(abs(diff_cpG_mapinfo)))) +
#  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
#  labs(title = "M0.1",x = "log10 (distance)",y = "Density") + theme_bw()+
#  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"),
#        axis.title.x = element_text(size = 8),  
#        axis.title.y = element_text(size = 8))+xlim(-2,4)+ylim(0,1)

hist_M0.1_pos_log_distance=ggplot(M0.1_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() , aes(x = log10(abs(diff_cpG_mapinfo)))) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+xlim(-2,4)+ylim(0,1)

hist_M0.2_pos_log_distance=ggplot(M0.2_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.1_pos_log_distance=ggplot(M1.1_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.2_pos_log_distance=ggplot(M1.2_pos%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M0.1_neg_log_distance=ggplot(M0.1_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M0.2_neg_log_distance=ggplot(M0.2_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

hist_M1.1_neg_log_distance=ggplot(M1.1_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (distance)", y = "Density") +theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+ xlim(-2,4)+ylim(0,1)

hist_M1.2_neg_log_distance=ggplot(M1.2_neg%>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs(), aes(x = log10(diff_cpG_mapinfo))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5))+xlim(-2,4)+ylim(0,1)

log_distance_layout_1=(hist_M0.1_pos_log_distance|hist_M1.1_pos_log_distance)/(hist_M0.2_pos_log_distance|hist_M1.2_pos_log_distance)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_layout_1

#ggsave("hist_log_distance_1.pdf")

log_distance_layout_2=(hist_M0.1_neg_log_distance|hist_M1.1_neg_log_distance)/(hist_M0.2_neg_log_distance|hist_M1.2_neg_log_distance)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_layout_2

#ggsave("hist_log_distance_2.pdf")

plot_grid(log_distance_layout_1,log_distance_layout_2,labels = "AUTO",label_size = 15, ncol = 1)
#ggsave(6*8inches) potrait
#plot_grid(log_distance_layout_1,log_distance_layout_2,labels = c("A (ER+)","B (ER-)"), ncol = 1)
##########


######################################################################################################################
### Distribution of the log 10 (distance) between the CpG in the methylation probe and the start position of a gene.
######################################################################################################################

hist_M0.1_pos_log_distance_start=ggplot(M0.1_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1", x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M0.2_pos_log_distance_start=ggplot(M0.2_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M1.1_pos_log_distance_start=ggplot(M1.1_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1", x = "log10 (distance)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M1.2_pos_log_distance_start=ggplot(M1.2_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (distance)",y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)


hist_M0.1_neg_log_distance_start=ggplot(M0.1_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (distance)", y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M0.2_neg_log_distance_start=ggplot(M0.2_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (distance)",y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M1.1_neg_log_distance_start=ggplot(M1.1_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1", x = "log10 (distance)",y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

hist_M1.2_neg_log_distance_start=ggplot(M1.2_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs(), aes(x = log10(abs(diff_mapinfo_geneStart)))) +
  geom_histogram(aes(y = ..density..), bins = 35, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2", x = "log10 (distance)",y = "Density")  + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+
  ylim(0,0.9)+xlim(-2.5,7.5)

log_distance_start_layout_1=(hist_M0.1_pos_log_distance_start|hist_M1.1_pos_log_distance_start)/(hist_M0.2_pos_log_distance_start|hist_M1.2_pos_log_distance_start)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_start_layout_1

#ggsave("hist_log_distance_start_1.pdf")

log_distance_start_layout_2=(hist_M0.1_neg_log_distance_start|hist_M1.1_neg_log_distance_start)/(hist_M0.2_neg_log_distance_start|hist_M1.2_neg_log_distance_start)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
log_distance_start_layout_2

#ggsave("hist_log_distance_start_2.pdf")

plot_grid(log_distance_start_layout_1,log_distance_start_layout_2,labels = "AUTO", label_size = 15,ncol = 1)
#ggsave(6*8inches)
############


######################################################################################################################
### Distribution of the log 10 (length) between the start position and the end position of a gene. 
######################################################################################################################

hist_M0.1_pos_logLength=ggplot(M0.1_pos%>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M0.2_pos_logLength=ggplot(M0.2_pos %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M1.1_pos_logLength=ggplot(M1.1_pos %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M1.2_pos_logLength=ggplot(M1.2_pos %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (length)",y = "Density")+ theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M0.1_neg_logLength=ggplot(M0.1_neg %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1", x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M0.2_neg_logLength=ggplot(M0.2_neg %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M1.1_neg_logLength=ggplot(M1.1_neg %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

hist_M1.2_neg_logLength=ggplot(M1.2_neg %>% select(gene_length) %>% distinct(), aes(x = log10(gene_length))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "log10 (length)",y = "Density") + theme_bw()+
  theme(plot.title = element_text(size=10,hjust = 0.5,face = "plain"))+ylim(0,0.8)+xlim(1,7)

logLength_layout_1=(hist_M0.1_pos_logLength|hist_M1.1_pos_logLength)/(hist_M0.2_pos_logLength|hist_M1.2_pos_logLength)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
logLength_layout_1

#ggsave("hist_logLength_1.pdf")

logLength_layout_2=(hist_M0.1_neg_logLength|hist_M1.1_neg_logLength)/(hist_M0.2_neg_logLength|hist_M1.2_neg_logLength)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
logLength_layout_2

#ggsave("hist_logLength_2.pdf")

plot_grid(logLength_layout_1,logLength_layout_2,labels = "AUTO", label_size = 15,ncol = 1)
#6*8 potrait
###########


######################################################################################################################
# Distribution of the relation to CpG island :five sections: Island, N shore, S shore, N shelf, S shelf, and No Info 
######################################################################################################################

## pie_M0.1_pos
pie_dat_M0.1_pos <- M0.1_pos %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M0.1_pos[1,1]="No Info"

## pie_M0.2_pos
pie_dat_M0.2_pos <- M0.2_pos %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M0.2_pos[1,1]="No Info"

## pie_M1.1_pos
pie_dat_M1.1_pos <- M1.1_pos %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M1.1_pos[1,1]="No Info"

## pie_M1.2_pos
pie_dat_M1.2_pos <- M1.2_pos %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M1.2_pos[1,1]="No Info"

pie_dat_M0.1_neg <- M0.1_neg %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M0.1_neg[1,1]="No Info"

pie_dat_M0.2_neg <- M0.2_neg %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M0.2_neg[1,1]="No Info"

pie_dat_M1.1_neg <- M1.1_neg %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M1.1_neg[1,1]="No Info"

pie_dat_M1.2_neg <- M1.2_neg %>% count(Relation_to_UCSC_CpG_Island) %>% mutate(percentage = n / sum(n) * 100)
pie_dat_M1.2_neg[1,1]="No Info"

#########funtion########
create_pie <- function(data_pie, title) {pie(data_pie$n, labels = paste0(data_pie$Relation_to_UCSC_CpG_Island," (",round(data_pie$percentage, 1), "%)"), col = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'), main = title,cex = 1,cex.main=1.5 )}

par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
#create_pie(pie_dat_M0.1_pos,"M0.1 pos")

pie_M0.1_pos <- create_pie(pie_dat_M0.1_pos,"M0.1")
pie_M1.1_pos <- create_pie(pie_dat_M1.1_pos,"M1.1")
pie_M0.2_pos <- create_pie(pie_dat_M0.2_pos,"M0.2")
pie_M1.2_pos <- create_pie(pie_dat_M1.2_pos,"M1.2")
mtext("ER+", outer = TRUE, cex = 1.4, font = 2)
mtext("A", side = 3, adj = 0, outer = TRUE, line = -1.5, cex = 2, font = 1)

#"pie_pos.pdf" 5*7 landscape
#pie_pos=(pie_M0.1_pos|pie_M1.1_pos)/(pie_M0.2_pos|pie_M1.2_pos)+
#  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
#pie_pos
#ggsave("pie_1_pos.pdf")

par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
create_pie(pie_dat_M0.1_neg,"M0.1")
create_pie(pie_dat_M1.1_neg,"M1.1")
create_pie(pie_dat_M0.2_neg,"M0.2")
create_pie(pie_dat_M1.2_neg,"M1.2")
mtext("ER-", outer = TRUE, cex = 1.4, font = 2)
mtext("B", side = 3, adj = 0, outer = TRUE, line = -1.5, cex = 2, font = 1)
#"pie_neg.pdf" 5*7 landscape
#ggsave("pie_2_neg.pdf")
############



######################################################################################################################
# Distribution of methylation mean values for the baseline and mediation models
######################################################################################################################

hist_M0.1_pos_MethyMean=ggplot(M0.1_pos, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M0.2_pos_MethyMean=ggplot(M0.2_pos, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M1.1_pos_MethyMean=ggplot(M1.1_pos, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M1.2_pos_MethyMean=ggplot(M1.2_pos, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M0.1_neg_MethyMean=ggplot(M0.1_neg, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M0.2_neg_MethyMean=ggplot(M0.2_neg, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M1.1_neg_MethyMean=ggplot(M1.1_neg, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "Mean methylation",y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

hist_M1.2_neg_MethyMean=ggplot(M1.2_neg, aes(x = Methyl_mean)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "Mean methylation", y = "Density") + theme_bw()+xlim(-5.0,5.0)+ylim(0,0.65)

MethyMean_layout_1=(hist_M0.1_pos_MethyMean|hist_M1.1_pos_MethyMean)/(hist_M0.2_pos_MethyMean|hist_M1.2_pos_MethyMean)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
MethyMean_layout_1

#ggsave("hist_MethylMean_1.pdf")

MethyMean_layout_2=(hist_M0.1_neg_MethyMean|hist_M1.1_neg_MethyMean)/(hist_M0.2_neg_MethyMean|hist_M1.2_neg_MethyMean)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
MethyMean_layout_2

#ggsave("hist_MethylMean_2.pdf")

plot_grid(MethyMean_layout_1,MethyMean_layout_2,labels = "AUTO",  label_size = 15,ncol = 1)
#6*8 potrait ("hist_MethylMean.pdf")
###########




######################################################################################################################
# Distribution of the GC content
######################################################################################################################
#The figure shows the distribution of the GC content which is the percentage of C and G among the four nucleotide base pairs (A, T, C, and G) in the baseline and mediation models across the cancer subtypes.


p_M0.1_pos=ggplot(M0.1_pos %>% select(Gene...GC.content), aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M0.2_pos=ggplot(M0.2_pos %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M1.1_pos=ggplot(M1.1_pos %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M1.2_pos=ggplot(M1.2_pos %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#f781bf", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M0.1_neg=ggplot(M0.1_neg %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.1",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M0.2_neg=ggplot(M0.2_neg %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M0.2",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.25)

p_M1.1_neg=ggplot(M1.1_neg %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.1",x = "GC content", y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.1)

p_M1.2_neg=ggplot(M1.2_neg %>% select(Gene...GC.content) , aes(x = Gene...GC.content)) +
  geom_histogram(aes(y = ..density..), binwidth = 2.5, fill = "#377eb8", color = "black", alpha = 0.7) +
  labs(title = "M1.2",x = "GC content",y = "Density") + theme_bw()+xlim(28,75)+ylim(0,0.25)

layout_1=(p_M0.1_pos|p_M1.1_pos)/(p_M0.2_pos|p_M1.2_pos)+
  plot_annotation(title = "ER+")& theme(plot.title = element_text(hjust = 0.5,size = 12))
layout_1
#ggsave("Figure_GC_content_1.pdf")

layout_2=(p_M0.1_neg|p_M1.1_neg)/(p_M0.2_neg|p_M1.2_neg)+
  plot_annotation(title = "ER-")& theme(plot.title = element_text(hjust = 0.5,size = 12))
layout_2
#ggsave("Figure_GC_content_2.pdf")

plot_grid(layout_1,layout_2,labels = "AUTO", label_size = 15,ncol = 1)
#6*8 potrait

#########



################### calculate average distance ############
# Mean of the log10 (distance) between the nearby CpG island and the CpG in the methylation probe.

M0.1_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M0.2_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M1.1_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M1.2_pos %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()

M0.1_neg %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M0.2_neg %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M1.1_neg %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()
M1.2_neg %>% select(diff_cpG_mapinfo) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_cpG_mapinfo > -2) %>% pull() %>% summary()

# Mean of the log 10 (distance) between the CpG in the methylation probe and the start position of a gene.
M0.1_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M0.2_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M1.1_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M1.2_pos %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()

M0.1_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M0.2_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M1.1_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()
M1.2_neg %>% select(diff_mapinfo_geneStart) %>% drop_na() %>% abs() %>% log10() %>% filter(diff_mapinfo_geneStart > -2) %>% pull() %>% summary()

# Mean of the log 10 (length) between the start position and the end position of a gene.
M0.1_pos %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M0.2_pos %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M1.1_pos %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M1.2_pos %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()

M0.1_neg %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M0.2_neg %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M1.1_neg %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()
M1.2_neg %>% select(gene_length) %>% drop_na() %>% abs() %>% log10() %>% filter(gene_length > -2) %>% pull() %>% summary()

# Average value of methylation mean values in the baseline and mediation models 
#M0.1_pos %>% select(Methyl_mean) %>% summary()

summary(M0.1_pos$Methyl_mean)
summary(M0.2_pos$Methyl_mean)
summary(M1.1_pos$Methyl_mean)
summary(M1.2_pos$Methyl_mean)

## GC content
summary(M0.1_pos$Gene...GC.content)
summary(M0.2_pos$Gene...GC.content)
summary(M1.1_pos$Gene...GC.content)
summary(M1.2_pos$Gene...GC.content)

summary(M0.1_neg$Gene...GC.content)
summary(M0.2_neg$Gene...GC.content)
summary(M1.1_neg$Gene...GC.content)
summary(M1.2_neg$Gene...GC.content)
