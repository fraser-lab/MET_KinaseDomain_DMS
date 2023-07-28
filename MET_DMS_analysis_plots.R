#####################################################################################################
# Code by: Gabriella Estevam @ UCSF
#
# analysis and plots of TPR-MET +/- Exon14 conditions
# Utilizes "dms_analysis_utilities.R" and "MET_V2.RMD" code and outputs written by Christian MacDonald @ UCSF 
#
# the goal of this code is to  analyze each experiment, compare the two conditions,
# and define statistical thresholds for comparison 
# this code uses the R bio3D package for structure mapping and computing (RMSD and Centroids)
#
# MET = MET delta Exon14 (Exon14 skipped ICD)
# Ex14 = WT MET (WT ICD)
#
#####################################################################################################


#####################################################################################################
###---------------------------------------- load libraries ------------------------------------------
#####################################################################################################

library(ggplot2)
library(gglorenz)
library(ggpubr)
library(readr)
library(stringr)
library(tidyverse)
library(dplyr)
library(ggridges)
library(ggbraid)
library(tidyquant)
library(bio3d)
source("dms_analysis_utilities.R")
library(scales)
library(colorspace)
library(ggpubr)
library(rstatix)
library(patchwork)
library(cowplot)
library(dendextend)
library("RColorBrewer")
library(gghighlight)
library(stringdist)

# color palette
cbp1 <- c("#999999","#D55E00", "#56B4E9","#009E73","#E69F00",
          "#0072B2","#F0E442","#D55E00", "#CC79A7")

#####################################################################################################
###------------------------------------ define initial data frames ----------------------------------
#####################################################################################################


require(data.table)
### writes this data frame into a permanent data file - this should be a ONE time run, then commented out

ex14_scores_ = ex14_scores
met_scores_ = met_scores

#write.csv(ex14_scores_, ".../ex14_scores.csv", row.names=FALSE)
#write.csv(met_scores, ".../met_scores.csv", row.names=FALSE)


### open data file as data frame for all raw inhibitor scores 
#ex14_scores_<- read.csv("ex14_scores.csv")
#met_scores_<- read.csv("met_scores.csv")

#met_pValues <-as.data.frame(fread("met_v2_main_identifiers_scores_pvalues_wt.tsv"))
#ex14_pValues <-as.data.frame(fread("ex14_v2_main_identifiers_scores_pvalues_wt.tsv"))


ex14_scores_ = ex14_scores_ %>% filter(mutation_type != "X")
met_scores_ = met_scores_ %>% filter(mutation_type != "X")

ex14_avg_scores <- ex14_scores_ %>% group_by(pos) %>% summarise(ex14_avg_IL3 = mean(IL3_score, na.rm=TRUE),
                                                               ex14_avg_IL3_withdrawal = mean(IL3_withdrawal_score, na.rm=TRUE),
                                                               ex14_avg_IL3_withdrawal_SE = mean(IL3_withdrawal_SE, na.rm=TRUE))

met_avg_scores <- met_scores_ %>% group_by(pos) %>% summarise(met_avg_IL3 = mean(IL3_score, na.rm=TRUE),
                                                             met_avg_IL3_withdrawal = mean(IL3_withdrawal_score, na.rm=TRUE),
                                                             met_avg_IL3_withdrawal_SE = mean(IL3_withdrawal_SE, na.rm=TRUE))

#### IL3 withdrawal 

ex_14_df<-data.frame(hgvs=ex14_scores_$hgvs,
                     ex14_SE = ex14_scores_$IL3_withdrawal_SE,
                     ex14_score=ex14_scores_$IL3_withdrawal_score,       
                     mutation_type=ex14_scores_$mutation_type,
                     pos=ex14_scores_$pos,
                     variants=ex14_scores_$variants)


met_df<-data.frame(hgvs=met_scores_$hgvs, 
                   met_SE = met_scores_$IL3_withdrawal_SE,
                   met_score =met_scores_$IL3_withdrawal_score,
                   mutation_type=met_scores_$mutation_type,
                   pos=met_scores_$pos, 
                   variants=met_scores_$variants)

#met_pValues_IL3_withdrawal_df <-data.frame(hgvs=met_pValues$hgvs,
#                                        met_IL3_withdrawal_pvalue = met_pValues$met_IL3_withdrawal_pvalue,
#                                        met_IL3_withdrawal_zscore = met_pValues$met_IL3_withdrawal_z) 

#ex14_pValues_IL3_withdrawal_df <-data.frame(hgvs=ex14_pValues$hgvs,
#                                         ex14_IL3_withdrawal_pvalue = ex14_pValues$ex14_IL3_withdrawal_pvalue,
#                                         ex14_IL3_withdrawal_zscore = ex14_pValues$ex14_IL3_withdrawal_z) 

#### IL3 

ex_14_df2<-data.frame(hgvs=ex14_scores_$hgvs,
                     ex14_SE = ex14_scores_$IL3_SE,
                     ex14_score_IL3 = ex14_scores_$IL3_score,
                     mutation_type_IL3 = ex14_scores_$mutation_type,
                     pos_IL3=ex14_scores_$pos,
                     variants=ex14_scores_$variants)

met_df2<-data.frame(hgvs=met_scores_$hgvs,
                   met_SE = met_scores_$IL3_SE,
                   met_score_IL3=met_scores_$IL3_score,
                   mutation_type_IL3=met_scores_$mutation_type,
                   pos_IL3=met_scores_$pos,
                   variants_IL3=met_scores_$variants)

#met_pValues_IL3_df <-data.frame(hgvs=met_pValues$hgvs,
#                                met_IL3_pvalue = met_pValues$met_IL3_pvalue,
#                                met_IL3_zscore = met_pValues$met_IL3_z)

#ex14_pValues_IL3_df <-data.frame(hgvs=ex14_pValues$hgvs,
#                                 ex14_IL3_pvalue = ex14_pValues$ex14_IL3_pvalue,
#                                 ex14_IL3_zscore = ex14_pValues$ex14_IL3_z) 


###---------------- merge data frames ---------------------------------

#met_df_pvalues <- merge(met_df, met_pValues_IL3_withdrawal_df)
#ex14_df_pvalues <- merge(ex_14_df, ex14_pValues_IL3_withdrawal_df)

#met_df2_pvalues <- merge(met_df2, met_pValues_IL3_df)
#ex14_df2_pvalues <- merge(ex_14_df2, ex14_pValues_IL3_df)

#merged_datasets_pvalue<-merge(ex14_df_pvalues,met_df_pvalues) ## IL3 withdrawal 
merged_datasets<-merge(ex_14_df,met_df)
merged_datasets_IL3<-merge(ex_14_df2,met_df2)

merged_datasets$score_diff <- abs((merged_datasets$met_score)-(merged_datasets$ex14_score)) #absolute difference 
merged_datasets$prop_error <- ((merged_datasets$met_SE)^2 +(merged_datasets$ex14_SE)^2)^(1/2) #propagated error
merged_datasets$raw_score_diff <- (merged_datasets$met_score)-(merged_datasets$ex14_score) #not absolute difference 
merged_datasets$ex14_abs_score <- abs(merged_datasets$ex14_score)
merged_datasets$met_abs_score <- abs(merged_datasets$met_score)

## adds propagation of error between scores / condition
## propagation of error is the square root of the sum of the squares of the errors in the quantities being added or subtracted

#merged_datasets_IL3_pvalue<-merge(ex14_df2_pvalues,met_df2_pvalues) ## IL3 
merged_datasets_avg<-merge(ex14_avg_scores,met_avg_scores )

## prepare data frames with statistically different and similar filters
scores_filtered_diff <- subset(merged_datasets, score_diff >= prop_error & mutation_type != "S"& mutation_type != "N") # statistically different 
scores_filtered_sim <- subset(merged_datasets, abs(score_diff) < 0.3 & abs(score_diff) >= prop_error & mutation_type != "S"& mutation_type != "N") #statistically similar 




#####################################################################################################
# --------------- Initial validation and characterization of mutational distributions ---------------
#####################################################################################################

###---------------- Plot IL3 mutations --------------------------------

# plot scatters of each mutation 
plot <- ggplot(merged_datasets_IL3, aes(x=met_score_IL3,y=ex14_score_IL3)) + 
  geom_point(shape=21, aes(colour = mutation_type_IL3))+
  theme_linedraw()+
  xlab("MET-Ex14") + ylab("MET+Ex14")
plot(plot)

plot <- ggplot(merged_datasets_IL3, aes(x=met_score_IL3, y=ex14_score_IL3))+
  geom_point(aes(colour=mutation_type_IL3))+
  #geom_density_2d()+
  theme_linedraw()+
  xlab("MET-Ex14") + ylab("MET+Ex14")
plot(plot)


# plot scatters with marginal distributions 
plot3 <- ggplot(merged_datasets_IL3 %>% filter(mutation_type_IL3 != "N"), aes(x=met_score_IL3, y=ex14_score_IL3, color = mutation_type_IL3)) + 
  geom_point(shape=1,aes(color = mutation_type_IL3)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens3 <- ggplot(merged_datasets_IL3 %>% filter(mutation_type_IL3 != "N"), aes(x = met_score_IL3, fill = mutation_type_IL3)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens4 <- ggplot(merged_datasets_IL3 %>% filter(mutation_type_IL3 != "N"), aes(x = ex14_score_IL3, fill = mutation_type_IL3)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))




###---------------- Plot IL3 withdrawal mutations -----------------------

# scatter
plot <- ggplot(data=merged_datasets, aes(x=met_score, y=ex14_score, colour = mutation_type)) + 
  geom_point(shape=21,aes(y=ex14_score, colour = mutation_type))+
  theme_linedraw()+
  xlab("MET-Ex14") + ylab("MET+Ex14")
plot(plot)

# scatter with marginal distribution
plot3 <- ggplot(merged_datasets %>% filter(mutation_type != "N"), aes(x=met_score, y=ex14_score, color = mutation_type, label= hgvs)) + 
  geom_point(shape=20,alpha=0.5,aes(color = mutation_type))+
  theme_pubr() +
  #geom_text()+
  theme(legend.position = "none")+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  geom_vline(xintercept=met_syn_mean, linetype="dashed")+
  geom_hline(yintercept=ex14_syn_mean,linetype="dashed")+
  geom_point(aes(x = -3.72715713, y =1.04762765), colour="black", size = 2)+  #S1122Q
  geom_point(aes(x = 0.861562751, y =-8.833887847), colour="black", size = 2)+ #L1062D 
  geom_point(aes(x =  -6.7885401, y =0.62601057), colour="black", size = 2)+ #V1121G
  geom_point(aes(x = -5.3728094, y =0.30282509), colour="black", size = 2)+ #V1121Y
  
  scale_color_manual(values=c("#219ebc", "#999999", "#56B4E9"))
dens3 <- ggplot(merged_datasets %>% filter(mutation_type != "N"), aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#219ebc", "#999999", "#56B4E9"))
dens4 <- ggplot(merged_datasets %>% filter(mutation_type != "N"), aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#219ebc", "#999999", "#56B4E9"))
dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))




#####################################################################################################
#--------------Distributions of syn, missense, and nonsense pre and post selection-------------------
#####################################################################################################

plot1 <- ggplot(ex14_scores %>% filter(mutation_type != "X" & mutation_type !="M" & mutation_type != "N")) +
  geom_histogram(aes(x = IL3_score), binwidth = 0.5, color = "#E69F00",fill="#E69F00")+
  xlab("Activity score") + 
  ylab("Density") + 
  xlim(-10,5)+
  theme_classic() + 
  geom_vline(xintercept = 0, linetype="dotted") +
  ggtitle("Synonymous IL3")+
  labs(color = "Mutation type")+
  theme(plot.title = element_text(size = 10, face = "bold"))


plot2 <- ggplot(ex14_scores %>% filter(mutation_type != "X" & mutation_type !="S" & mutation_type != "N")) +
  geom_histogram(aes(x = IL3_score), binwidth = 0.5, color = "#56B4E9",fill="#56B4E9") +
  xlab("Activity score") + 
  ylab("Density") + 
  xlim(-10,5)+
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic() + 
  ggtitle("Missense IL3")+
  labs(color = "Mutation type")+
  theme(plot.title = element_text(size = 10, face = "bold"))


plot3 <- ggplot(ex14_scores %>% filter(mutation_type != "X" & mutation_type !="S" & mutation_type != "M")) +
  geom_histogram(aes(x = IL3_score), binwidth = 0.5, color = "#009E73",fill="#009E73") +
  xlab("Activity score") + 
  ylab("Density") + 
  xlim(-10,5)+
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic() + 
  ggtitle("Nonsense IL3")+
  labs(color = "Mutation type")+
  theme(plot.title = element_text(size = 10, face = "bold"))



plot4 <- ggplot(ex14_scores %>% filter(mutation_type != "X" & mutation_type !="M" & mutation_type != "N")) +
  geom_histogram(aes(x = IL3_withdrawal_score), binwidth = 0.5,color = "#E69F00",fill="#E69F00") +
  xlab("Activity score") + 
  ylab("Density") + 
  xlim(-10,5)+
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic() + 
  ggtitle("Synonymous IL3 Withdrawal")+
  labs(color = "Mutation type")+
  theme(plot.title = element_text(size = 10, face = "bold"))

plot5 <- ggplot(ex14_scores %>% filter(mutation_type != "X" & mutation_type != "S" & mutation_type != "N")) +
  geom_histogram(aes(x = IL3_withdrawal_score), binwidth = 0.5,color = "#56B4E9",fill="#56B4E9" ) +
  xlab("Activity score") + 
  ylab("Density") + 
  xlim(-10,5)+
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic() + 
  ggtitle("Missense IL3 Withdrawal")+
  labs(color = "Mutation type")+
  theme(plot.title = element_text(size = 10, face = "bold"))

plot6 <- ggplot(ex14_scores %>% filter(mutation_type != "X" & mutation_type != "M" & mutation_type != "S")) +
  geom_histogram(aes(x = IL3_withdrawal_score), bins = 5, binwidth = 0.5,color = "#009E73",fill="#009E73" ) +
  xlab("Activity score") + 
  ylab("Density") + 
  xlim(-10,5)+
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic() + 
  ggtitle("Nonsense IL3 Withdrawal")+
  labs(color = "Mutation type")+
  theme(plot.title = element_text(size = 10, face = "bold"))

plot_grid(plot1,plot4,plot2,plot5,plot3,plot6, ncol=2, nrow=3)



#####################################################################################################
#-----------------------Comparison Plots of statistically different mutations -----------------------
#####################################################################################################

# data frame of scores that are statistically different 
met_sum_scores <- scores_filtered_diff %>% group_by(pos) %>% summarise(met= sum(met_score, na.rm=TRUE))
ex14_sum_scores <- scores_filtered_diff %>% group_by(pos) %>% summarise(ex14= sum(ex14_score, na.rm=TRUE))

# plot sums of scored that are statistically different
met_sum_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_sum_3R7O, met_sum_scores ,"met")
write.pdb(x, file="met_sum_3R7O.pdb")
ex14_sum_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_sum_3R7O, ex14_sum_scores ,"met")
write.pdb(x, file="ex14_sum_3R7O.pdb")

# plot statistically relevant DIFFERENT mutations
plot1 <- ggplot(scores_filtered_diff %>% filter(mutation_type != "N")) + 
  geom_point(shape=1, aes(x=met_score,y=ex14_score, color = mutation_type)) +
  #stat_smooth(method = "lm", fullrange = TRUE,se = FALSE) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens1 <- ggplot(scores_filtered_diff%>% filter(mutation_type != "N"), aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens2 <- ggplot(scores_filtered_diff%>% filter(mutation_type != "N"), aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens1 + plot_spacer() + plot1 + dens2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

dens1 + plot_spacer() + plot1 + dens2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))


####################################################################################################
####----------------------------------- proline 1153 analysis -------------------------------------
####################################################################################################

proline_site_distributions <- scores_filtered_diff %>% filter((pos %in% c(1153,1121,1124)) & mutation_type!= "S"& mutation_type!= "N")

plot_pro_site_diff <- ggplot(proline_site_distributions) + 
  geom_point(shape=1, aes(x=met_score,y=ex14_score, color = mutation_type, label=pos)) +
  #stat_smooth(method = "lm", fullrange = TRUE,se = FALSE) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  geom_vline(xintercept = 0, linetype = "dotted", color = "red")+
  geom_text(aes(x=met_score,y=ex14_score,label=pos),hjust=-.1, vjust=0)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
plot(plot_pro_site_diff )

P1153 = merged_datasets %>% filter(mutation_type != "N" & mutation_type != "S" & (pos %in% c(1153)))

plot_P1153<- ggplot()+
  geom_histogram((aes(x=P1153$ex14_score)),color="black",fill="grey", alpha=0.5)+
  geom_histogram((aes(x=P1153$met_score)),color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, color = "grey", linetype="dashed")+
  geom_vline(xintercept = met_syn_mean, color="deepskyblue2", linetype="dashed")+
  ggtitle("P1153 B5-turn site")+
  theme_bw() +
  theme(
    panel.background=element_rect(colour="black",size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank())
plot(plot_P1153 )


hydrophobic<- c('I','L','V','M','A','H','K')
positive<- c('R','K','H')
negative<- c('D','E')
polar_uncharged<-c('N','Q','T','S','C')
aromatic<-c('W','Y','F')

P1153_physio<-bind_rows(
  {P1153 %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {P1153 %>% filter(variants %in% negative)%>% group_by(variants, ID='negative')},
  {P1153 %>% filter(variants %in% positive)%>% group_by(variants, ID='positive')},
  {P1153 %>% filter(variants %in% polar_uncharged)%>% group_by(variants, ID='polar uncharged')},
  {P1153 %>% filter(variants %in% aromatic)%>% group_by(variants, ID='aromatic')}
)

plot_P1153_physio <- ggplot(P1153_physio, aes(x=factor(pos), y=ex14_score, fill = ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.6) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("")+
  theme_bw()+
  theme(text = element_text(size = 20),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_P1153_physio)

#####################################################################################################
#-----------------------Comparison Plots of statistically similar mutations -------------------------
#####################################################################################################

# datafrme of scores that are statistically SIMILAR
met_sum_sim_scores <- scores_filtered_sim %>% filter (mutation_type != "S") %>% group_by(pos) %>% summarise(met=abs(sum(met_score, na.rm=TRUE)))
ex14_sum_sim_scores <- scores_filtered_sim  %>% filter (mutation_type != "S")%>% group_by(pos) %>% summarise(ex14=abs(sum(ex14_score, na.rm=TRUE)))


#white the avg between +/-Ex14, since we already know the difference is small 

# plot sums of scored that are statistically SIMILAR
met_sum_sim_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_sum_sim_3R7O, met_sum_sim_scores ,"met")
write.pdb(x, file="met_sum_sim_3R7O.pdb")
ex14_sum_sim_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_sum_sim_3R7O, ex14_sum_sim_scores ,"ex14")
write.pdb(x, file="ex14_sum_sim_3R7O.pdb")

# plot statistically relevant SIMILAR mutations
plot1 <- ggplot(scores_filtered_sim %>% filter(mutation_type != "N" & mutation_type != "S")) + 
  geom_point(shape=1, aes(x=met_score,y=ex14_score, color = mutation_type)) +
  #stat_smooth(method = "lm", fullrange = TRUE,se = FALSE) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens1 <- ggplot(scores_filtered_sim%>% filter(mutation_type != "N" & mutation_type != "S"), aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens2 <- ggplot(scores_filtered_sim%>% filter(mutation_type != "N" & mutation_type != "S"), aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dens1 + plot_spacer() + plot1 + dens2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

#####################################################################################################
###------------------------------------ Replicate validation  ----------------------------------
#####################################################################################################



## MET replicates +IL3
ex14_rep1_IL3<-as.data.frame(fread(".../Ex14/tsv/R1_sel/main_identifiers_scores.tsv"))
ex14_rep2_IL3<-as.data.frame(fread(".../Ex14/tsv/R2_sel/main_identifiers_scores.tsv"))
ex14_rep3_IL3<-as.data.frame(fread(".../Ex14/tsv/R3_sel/main_identifiers_scores.tsv"))

ex14_hgvs_pos <- data.frame(hgvs=ex14_scores_$hgvs, pos = ex14_scores_$pos)

ex14_R1_IL3 <- data.frame(hgvs=ex14_rep1_IL3$hgvs, R1=ex14_rep1_IL3$score)
ex14_R1_IL3  = merge(ex14_R1_IL3,ex14_hgvs_pos)
ex14_R2_IL3 <- data.frame(hgvs=ex14_rep2_IL3$hgvs, R2=ex14_rep2_IL3$score)
ex14_R2_IL3  = merge(ex14_R2_IL3,ex14_hgvs_pos)
ex14_R3_IL3 <- data.frame(hgvs=ex14_rep3_IL3$hgvs, R3=ex14_rep3_IL3$score)
ex14_R3_IL3  = merge(ex14_R3_IL3,ex14_hgvs_pos)

ex14_R1_IL3_avg <- ex14_R1_IL3 %>% group_by(pos) %>% summarise(R1=mean(R1))
ex14_R2_IL3_avg <- ex14_R2_IL3 %>% group_by(pos) %>% summarise(R2=mean(R2))
ex14_R3_IL3_avg <- ex14_R3_IL3 %>% group_by(pos) %>% summarise(R3=mean(R3))

ex14_rep12_df_IL3 <- merge(ex14_R1_IL3_avg,ex14_R2_IL3_avg)
ex14_rep13_df_IL3 <- merge(ex14_R1_IL3_avg, ex14_R3_IL3_avg)

corr_ex14_rep12_df_IL3 <- cor(ex14_rep12_df_IL3$R1, ex14_rep12_df_IL3$R2, method = c("pearson"))
corr_ex14_rep13_df_IL3 <- cor(ex14_rep13_df_IL3$R1, ex14_rep13_df_IL3$R3, method = c("pearson"))


## MET replicates IL3 withdrawal
ex14_rep1_<-as.data.frame(fread(".../Ex14/tsv/W_R1_sel/main_identifiers_scores.tsv"))
ex14_rep2_<-as.data.frame(fread(".../Ex14/tsv/W_R2_sel/main_identifiers_scores.tsv"))
ex14_rep3_<-as.data.frame(fread(".../Ex14/tsv/W_R3_sel/main_identifiers_scores.tsv"))

ex14_R1 <- data.frame(hgvs=ex14_rep1_$hgvs, R1=ex14_rep1_$score)
ex14_R1 = merge(ex14_R1,ex14_hgvs_pos)
ex14_R2 <- data.frame(hgvs=ex14_rep2_$hgvs, R2=ex14_rep2_$score)
ex14_R2 = merge(ex14_R2,ex14_hgvs_pos)
ex14_R3 <- data.frame(hgvs=ex14_rep3_$hgvs, R3=ex14_rep3_$score)
ex14_R3 = merge(ex14_R3,ex14_hgvs_pos)

ex14_R1_avg <- ex14_R1 %>% group_by(pos) %>% summarise(R1=mean(R1))
ex14_R2_avg <- ex14_R2 %>% group_by(pos) %>% summarise(R2=mean(R2))
ex14_R3_avg <- ex14_R3 %>% group_by(pos) %>% summarise(R3=mean(R3))

ex14_rep12_df <- merge(ex14_R1_avg,ex14_R2_avg)
ex14_rep13_df <- merge(ex14_R1_avg,ex14_R3_avg)

corr_ex14_rep12_df <- cor(ex14_rep12_df$R1, ex14_rep12_df$R2, method = c("pearson"))
corr_ex14_rep13_df <- cor(ex14_rep13_df$R1, ex14_rep13_df$R3, method = c("pearson"))

plot_Ex14_R1_R2_IL3 <- ggplot(data = ex14_rep12_df_IL3, aes(x=R1, y=R2)) +
  ggtitle("(+IL3) MET R1 vs R2")+
  xlab("Activity score R1")+
  ylab("Activity score R2")+
  xlim(-6,6)+
  ylim(-6,6)+
  theme_classic()+
  geom_point()
plot_Ex14_R1_R3_IL3 <- ggplot(data = ex14_rep13_df_IL3, aes(x=R1, y=R3)) +
  ggtitle("(+IL3) MET R1vs R3")+xlab("Activity score R1")+
  ylab("Activity score R3")+
  xlim(-6,6)+
  ylim(-6,6)+
  theme_classic()+
  geom_point()+
  theme_classic()
plot_grid(plot_Ex14_R1_R2_IL3,plot_Ex14_R1_R3_IL3)

plot_Ex14_R1_R2 <- ggplot(data = ex14_rep12_df, aes(x=R1, y=R2)) +
  ggtitle("(-IL3) MET R1 vs R2")+
  xlab("Activity score R1")+
  ylab("Activity score R2")+
  theme_classic()+
  geom_point()
plot_Ex14_R1_R3 <- ggplot(data = ex14_rep13_df, aes(x=R1, y=R3)) +
  ggtitle("(-IL3) MET R1vs R3")+xlab("Activity score R1")+
  ylab("Activity score R3")+
  geom_point()+
  theme_classic()

plot_grid(plot_Ex14_R1_R2_IL3,plot_Ex14_R1_R3_IL3,plot_Ex14_R1_R2,plot_Ex14_R1_R3)



#### METdEx14 replicates +IL3
## MET replicates +IL3
met_rep1_IL3<-as.data.frame(fread(".../Met/tsv/R1_sel/main_identifiers_scores.tsv"))
met_rep2_IL3<-as.data.frame(fread(".../Met/tsv/R2_sel/main_identifiers_scores.tsv"))
met_rep3_IL3<-as.data.frame(fread(".../Met/tsv/R3_sel/main_identifiers_scores.tsv"))

met_hgvs_pos <- data.frame(hgvs=met_scores_$hgvs, pos = met_scores_$pos)

met_R1_IL3 <- data.frame(hgvs=met_rep1_IL3$hgvs, R1=met_rep1_IL3$score)
met_R1_IL3  = merge(met_R1_IL3,met_hgvs_pos)
met_R2_IL3 <- data.frame(hgvs=met_rep2_IL3$hgvs, R2=met_rep2_IL3$score)
met_R2_IL3  = merge(met_R2_IL3,met_hgvs_pos)
met_R3_IL3 <- data.frame(hgvs=met_rep3_IL3$hgvs, R3=met_rep3_IL3$score)
met_R3_IL3  = merge(met_R3_IL3,met_hgvs_pos)

met_R1_IL3_avg <- met_R1_IL3 %>% group_by(pos) %>% summarise(R1=mean(R1))
met_R2_IL3_avg <- met_R2_IL3 %>% group_by(pos) %>% summarise(R2=mean(R2))
met_R3_IL3_avg <- met_R3_IL3 %>% group_by(pos) %>% summarise(R3=mean(R3))

met_rep12_df_IL3 <- merge(met_R1_IL3_avg,met_R2_IL3_avg)
met_rep13_df_IL3 <- merge(met_R1_IL3_avg, met_R3_IL3_avg)

corr_met_rep12_df_IL3 <- cor(met_rep12_df_IL3$R1, met_rep12_df_IL3$R2, method = c("pearson"))
corr_met_rep13_df_IL3 <- cor(met_rep13_df_IL3$R1, met_rep13_df_IL3$R3, method = c("pearson"))


## MET replicates IL3 withdrawal
met_rep1_<-as.data.frame(fread(".../Met/tsv/W_R1_sel/main_identifiers_scores.tsv"))
met_rep2_<-as.data.frame(fread(".../Met/tsv/W_R2_sel/main_identifiers_scores.tsv"))
met_rep3_<-as.data.frame(fread(".../Met/tsv/W_R3_sel/main_identifiers_scores.tsv"))

met_R1 <- data.frame(hgvs=met_rep1_$hgvs, R1=met_rep1_$score)
met_R1 = merge(met_R1,met_hgvs_pos)
met_R2 <- data.frame(hgvs=met_rep2_$hgvs, R2=met_rep2_$score)
met_R2 = merge(met_R2,met_hgvs_pos)
met_R3 <- data.frame(hgvs=met_rep3_$hgvs, R3=met_rep3_$score)
met_R3 = merge(met_R3,met_hgvs_pos)

met_R1_avg <- met_R1 %>% group_by(pos) %>% summarise(R1=mean(R1))
met_R2_avg <- met_R2 %>% group_by(pos) %>% summarise(R2=mean(R2))
met_R3_avg <- met_R3 %>% group_by(pos) %>% summarise(R3=mean(R3))

met_rep12_df <- merge(met_R1_avg,met_R2_avg)
met_rep13_df <- merge(met_R1_avg,met_R3_avg)

corr_met_rep12_df <- cor(met_rep12_df$R1, met_rep12_df$R2, method = c("pearson"))
corr_met_rep13_df <- cor(met_rep13_df$R1, met_rep13_df$R3, method = c("pearson"))

plot_met_R1_R2_IL3 <- ggplot(data = met_rep12_df_IL3, aes(x=R1, y=R2)) +
  ggtitle("(+IL3) METdEx14 R1 vs R2")+
  xlab("Activity score R1")+
  ylab("Activity score R2")+
  xlim(-8,8)+
  ylim(-8,8)+
  theme_classic()+
  geom_point()
plot_met_R1_R3_IL3 <- ggplot(data = met_rep13_df_IL3, aes(x=R1, y=R3)) +
  ggtitle("(+IL3) METdEx14 R1vs R3")+xlab("Activity score R1")+
  ylab("Activity score R3")+
  xlim(-8,8)+
  ylim(-8,8)+
  theme_classic()+
  geom_point()+
  theme_classic()
plot_grid(plot_met_R1_R2_IL3,plot_met_R1_R3_IL3)

plot_met_R1_R2 <- ggplot(data = met_rep12_df, aes(x=R1, y=R2)) +
  ggtitle("(-IL3) METdEx14 R1 vs R2")+
  xlab("Activity score R1")+
  ylab("Activity score R2")+
  theme_classic()+
  geom_point()
plot_met_R1_R3 <- ggplot(data = met_rep13_df, aes(x=R1, y=R3)) +
  ggtitle("(-IL3) METdEx14 R1 vs R3")+xlab("Activity score R1")+
  ylab("Activity score R3")+
  geom_point()+
  theme_classic()

plot_grid(plot_met_R1_R2_IL3,plot_met_R1_R3_IL3,plot_met_R1_R2,plot_met_R1_R3)


#####################################################################################################
###--------------------------- plot syn mutational distributions-------------------------------------
#####################################################################################################

met_syn = (met_scores_ %>% filter(mutation_type != "M" & mutation_type != "N" &mutation_type != "X"))
ex14_syn = (ex14_scores_ %>% filter(mutation_type != "M" & mutation_type != "N" &mutation_type != "X"))

plot_met_wt_dist<-ggplot(met_syn, aes(x = IL3_withdrawal_score)) + 
  geom_histogram(aes(y =..density..), 
                 colour = "black", 
                 fill = "#56B4E9",
                 alpha=.2) +
  ggtitle("METdEx14 WT ditribution")+
  geom_vline(xintercept=mean(met_syn$IL3_withdrawal_score), size=1, color="red", linetype="dotted")+
  xlim(-5,5)

plot_ex14_wt_dist<-ggplot(ex14_syn, aes(x = IL3_withdrawal_score)) + 
  geom_histogram(aes(y =..density..), 
                 colour = "black", 
                 fill = "#56B4E9",
                 alpha=.2) +
  ggtitle("MET+Ex14 WT ditribution")+
  geom_vline(xintercept=mean(ex14_syn$IL3_withdrawal_score), size=1, color="red", linetype="dotted")+
  xlim(-5,5)

plot_grid(plot_ex14_wt_dist,plot_met_wt_dist,ncol=2,nrow=1)


met_syn_mean = mean(met_syn$IL3_withdrawal_score) # -0.4341432
ex14_syn_mean = mean(ex14_syn$IL3_withdrawal_score) # -0.1426368


#####################################################################################################
#----------------------- statistical cut offs for GOF, LOF mutations  -------------------------------
#####################################################################################################

# get the spread of spread the data and WT syn 
# get spread of standard error of WT and NS 
# if SE of mutant is <= mean SE of WT && the mutant is in the +/-95%-tile (or +/-2 SD from mean of WT) then use

## means of the SE for each condition
met_wt_SE_mean = mean(met_syn$IL3_withdrawal_SE) # MET_avg_SE = 0.4603832
ex14_wt_SE_mean = mean(ex14_syn$IL3_withdrawal_SE) # Ex14_avg_SE = 0.2993557

### plots the WT SE distributions 
# METdEx14
plot_met_wt_SE<-ggplot(met_syn, aes(x = IL3_withdrawal_SE)) + 
  geom_histogram(aes(y =..density..), colour = "black",  fill = "#56B4E9", alpha=.2) +
  ggtitle("METdEx14 WT SE ditribution")+
  geom_vline(xintercept=mean(met_syn$IL3_withdrawal_SE), size=1, color="red", linetype="dotted")+
  xlim(-5,5)
plot(plot_met_wt_SE)

# MET+EX14
plot_ex14_wt_SE<-ggplot(ex14_syn, aes(x = IL3_withdrawal_SE)) + 
  geom_histogram(aes(y =..density..), colour = "black",  fill = "#56B4E9", alpha=.2) +
  ggtitle("MET+Ex14 WT SE distribution")+
  geom_vline(xintercept=mean(ex14_syn$IL3_withdrawal_SE), size=1, color="red", linetype="dotted")+
  xlim(-5,5)
plot(plot_ex14_wt_SE)

plot_grid(plot_ex14_wt_SE,plot_met_wt_SE,ncol=2,nrow=1)

#### some new definitions
met_missense = (met_scores_ %>% filter(mutation_type != "X" & mutation_type !="S" & mutation_type != "N"))
ex14_missense = (ex14_scores_ %>% filter(mutation_type != "X" & mutation_type !="S" & mutation_type != "N"))
met_2SD <- subset(met_missense%>% filter(abs(IL3_withdrawal_score - mean(met_syn$IL3_withdrawal_score)) >= 2*sd(met_syn$IL3_withdrawal_score)))
ex14_2SD <- subset(ex14_missense%>% filter(abs(IL3_withdrawal_score - mean(ex14_syn$IL3_withdrawal_score)) >= 2*sd(ex14_syn$IL3_withdrawal_score)))

### plots the mutations that are +/- 2SD 
plot_met_missense <- ggplot(met_missense) +
  geom_histogram(aes(x = IL3_withdrawal_score)) +
  xlab("Activity score") + 
  ylab("Count") + 
  xlim(-10,10)+
  geom_vline(xintercept = -0.4341432, linetype="dotted", color = "red") +
  ggtitle("Missense Mutation Distribution")+
  
plot(plot_met_missense)

plot_Overlay_distributions_missense <- ggplot() +
  geom_histogram(aes(x=met_missense$IL3_withdrawal_score), color = "deepskyblue", fill = "deepskyblue", alpha = 0.3) + 
  geom_vline(xintercept = 2*sd(met_syn$IL3_withdrawal_score),linetype="dashed", color = "skyblue", size=1) +
  theme_classic()+
  ggtitle("METdEx14")+
  xlab("Activity Score")+
  ylim(0,1000)+
  ylab("Count")+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_Overlay_distributions_missense)

plot_Overlay_distributions_missense2 <- ggplot() +
  geom_histogram(aes(x=ex14_missense$IL3_withdrawal_score), color = "orange", fill = "orange", alpha = 0.2) + 
  geom_vline(xintercept = 2*sd(ex14_syn$IL3_withdrawal_score),linetype="dashed", color = "orange", size=1) +
  theme_classic()+
  ylim(0,1000)+
  ggtitle("MET")+
  xlab("Activity Score")+
  ylab("Count")+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_Overlay_distributions_missense)

plot_grid(plot_Overlay_distributions_missense,plot_Overlay_distributions_missense2, ncol=2, nrow=1)

#####################################################################################################
#---------------------- Histograms of of GOF METdEx14 and LOF +Ex14 ---------------------------------
#####################################################################################################

#these are comparative histograms showing statistically significant score differences based on SE progigation of error 
# takes into actount: 
# 1) propagation of error 
# 2) mutations that are 2 standard deviations away from the mean of their respective syn distribution 

met_2SD_df <- data.frame(hgvs=met_2SD$hgvs,
                         met_SE = met_2SD$IL3_withdrawal_SE,
                         met_score=met_2SD$IL3_withdrawal_score,       
                         mutation_type=met_2SD$mutation_type,
                         pos=met_2SD$pos,
                         variants=met_2SD$variants)

ex14_2SD_df <- data.frame(hgvs=ex14_2SD$hgvs,
                         ex14_SE = ex14_2SD$IL3_withdrawal_SE,
                         ex14_score=ex14_2SD$IL3_withdrawal_score,       
                         mutation_type=ex14_2SD$mutation_type,
                         pos=ex14_2SD$pos,
                         variants=ex14_2SD$variants)

merge1 <- merge(scores_filtered_diff, met_2SD_df)
merge2 <- merge(scores_filtered_diff, ex14_2SD_df)
merge3 <- merge(merge1,merge2)

met_GOF_dist <- subset(merge3 %>% filter(met_score>=0))
ex14_LOF_dist <- subset(merge3 %>% filter(ex14_score<=0))

met_LOF_dist <- subset(merge3 %>% filter(met_score<=0))
ex14_GOF_dist <- subset(merge3 %>% filter(ex14_score>=0))

# creates data frame of positions that have the most divergent activity score differences
metGOF_ex14LOF_same_pos <- merge(met_GOF_dist,ex14_LOF_dist)
metLOF_ex14GOF_same_pos <- merge(met_LOF_dist,ex14_GOF_dist)

# plots scatter of METdEx14 GOF
met_2SD_df_filtered = met_2SD_df %>% filter (met_score>0)
ex14_2SD_df_filtered = ex14_2SD_df %>% filter (ex14_score>0)

plot_met_GOF_scatter <- ggplot() +
  geom_point(aes(x=met_2SD_df_filtered$pos, y=met_2SD_df_filtered$met_score), size=3, alpha=0.7, color="deepskyblue4")+
  xlab("Position")+
  ylab("Activity Score")+
  ggtitle("METdEx14 GOF")+
  xlim(1059,1345)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_met_GOF_scatter)

plot_ex14_GOF_scatter <- ggplot() +
  geom_point(aes(x=ex14_2SD_df_filtered$pos, y=ex14_2SD_df_filtered$ex14_score), size=3, alpha=0.7, color="deepskyblue4")+
  xlab("Position")+
  ylab("Activity Score")+
  ggtitle("MET GOF")+
  xlim(1059,1345)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_ex14_GOF_scatter)
plot_grid(plot_met_GOF_scatter,plot_ex14_GOF_scatter, ncol = 2, nrow=1)

### this is what is in the paper figure 5
#met_GOF_avg_scores <- met_2SD_df_filtered %>% group_by(pos) %>% summarise(met_GOF_avg = mean(met_score, na.rm=TRUE))
met_GOF_avg_scores <- met_2SD_df_filtered %>% group_by(pos) %>% summarise(met_GOF_count = n()) 
met_GOF_avg_scores_3R7O <- read.pdb("3R7O")
#x = map_scores_pdb(met_GOF_avg_scores_3R7O, met_GOF_avg_scores, "met_GOF_avg")
x = map_scores_pdb(met_GOF_avg_scores_3R7O, met_GOF_avg_scores, "met_GOF_count")
#write.pdb(x, file="met_GOF_avg_scores_3R7O.pdb")
write.pdb(x, file="met_GOF_count_3R7O.pdb")

#ex14_GOF_avg_scores <- ex14_2SD_df_filtered %>% group_by(pos) %>% summarise(ex14_GOF_avg = mean(ex14_score, na.rm=TRUE))
ex14_GOF_avg_scores <- ex14_2SD_df_filtered %>% group_by(pos) %>% summarise(ex14_GOF_count = n())
ex14_GOF_avg_scores_3R7O <- read.pdb("3R7O")
#x = map_scores_pdb(ex14_GOF_avg_scores_3R7O, ex14_GOF_avg_scores, "ex14_GOF_avg")
x = map_scores_pdb(ex14_GOF_avg_scores_3R7O, ex14_GOF_avg_scores, "ex14_GOF_count")
#write.pdb(x, file="ex14_GOF_avg_scores_3R7O.pdb")
write.pdb(x, file="ex14_GOF_count_3R7O.pdb")

#### GOF and LOF for figure 5


met_GOF_LOF_avg_scores <- met_GOF_LOF_avg_scores_ %>% 
  dplyr::group_by(pos) %>% 
  dplyr::summarise(met_GOF_LOF_avg = mean(met_score), counts = n()) %>% dplyr::ungroup() 

met_GOF_LOF_avg_scores_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_GOF_LOF_avg_scores_3R7O, met_GOF_LOF_avg_scores, "met_GOF_LOF_avg")
y = map_counts_pdb(met_GOF_LOF_avg_scores_3R7O, met_GOF_LOF_avg_scores, "counts")
write.pdb(x,y, file="met_GOF_LOF_avg_scores_3R7O.pdb")

met_GOF_LOF_avg_scores_3R7O_counts <- read.pdb("met_GOF_LOF_avg_scores_3R7O.pdb")
x = map_scores_pdb(met_GOF_LOF_avg_scores_3R7O_counts, met_GOF_LOF_avg_scores, "counts")
write.pdb(x, file="met_GOF_LOF_avg_scores_3R7O_counts.pdb")




ex14_GOF_LOF_avg_scores <- ex14_2SD_df %>% dplyr::group_by(pos) %>% dplyr::summarise(ex14_GOF_LOF_avg = mean(ex14_score, na.rm=TRUE))%>% dplyr::ungroup()
ex14_GOF_LOF_avg_scores_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_GOF_LOF_avg_scores_3R7O , ex14_GOF_LOF_avg_scores, "ex14_GOF_LOF_avg")
write.pdb(x, file="ex14_GOF_LOF_avg_scores_3R7O.pdb")



# GOF for both 
merged_GOF_for_both <- merge(ex14_2SD_df_filtered, met_2SD_df_filtered) 

#### plots histograms of positional GOF LOF 
plot_scores_diff3 <- ggplot() +
  geom_histogram(aes(x=metGOF_ex14LOF_same_pos$met_score),color="black", fill="deepskyblue4", alpha=0.6)+
  geom_histogram(aes(x=metGOF_ex14LOF_same_pos$ex14_score),color="black", fill="grey", alpha=0.6)+
  theme_classic()+
  xlab("Activity Score")+
  ggtitle("GOF METdEx14, LOF MET+Ex14")+
  xlim(-10,3)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_scores_diff4 <- ggplot() +
  geom_histogram(aes(x=metLOF_ex14GOF_same_pos$met_score),color="black", fill="deepskyblue4", alpha=0.6)+
  geom_histogram(aes(x=metLOF_ex14GOF_same_pos$ex14_score),color="black", fill="grey", alpha=0.6)+
  theme_classic()+
  xlab("Activity Score")+
  ggtitle("LOF METdEx14, GOF MET+Ex14")+ theme_bw()+
  xlim(-10,3)+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_grid(plot_scores_diff3,plot_scores_diff4,ncol=1,nrow=2)        

###### plots scatters of positional GOF LOF 

plot_scores_diff5 <- ggplot() +
  geom_point(aes(x=metGOF_ex14LOF_same_pos$met_score, y=metGOF_ex14LOF_same_pos$ex14_score), size=3, alpha=0.7, color="deepskyblue4")+
  xlab("METdEx14 Activity Score")+
  ylab("MET+Ex14 Activity Score")+
  ggtitle("GOF METdEx14, LOF MET+Ex14")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_scores_diff6 <- ggplot() +
  geom_point(aes(x=metLOF_ex14GOF_same_pos$met_score, y=metLOF_ex14GOF_same_pos$ex14_score),size=3, alpha=0.7,color="deepskyblue4")+
  ylab("MET+Ex14 Activity Score")+
  ggtitle("LOF METdEx14, GOF MET+Ex14")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_grid(plot_scores_diff5,plot_scores_diff6,ncol=1,nrow=2)   

#####################################################################################################
#------------------------------- Largest difference line plot --------------------------------------
#####################################################################################################

# scores_filtered_diff = data frame, already filtered by propigation of error
# score_diff = y-axis 
# position = x-axis

#avg_score_diff = mean(scores_filtered_diff$score_diff)
#avg_score_diff_SD = 4*sd(scores_filtered_diff$score_diff)
raw_scores_filtered_diff <- subset(merged_datasets, score_diff >= prop_efrror & mutation_type != "S"& mutation_type != "N") # statistically different 

## prepare data frames with statistically different and similar filters
avg_diff = raw_scores_filtered_diff %>% group_by(pos) %>% summarise(mean_diff = mean(raw_score_diff))
avg_score_diff = mean(avg_diff$mean_diff)
avg_score_diff_SD =sd(avg_diff$mean_diff)

avg_score_MOST_diff = avg_diff  %>% filter (mean_diff >= (avg_score_diff+avg_score_diff_SD))


plot_met_ex_pos_difference <- ggplot()+
  geom_line(aes(x=avg_diff$pos, y=avg_diff$mean_diff), color = "#219ebc", size=0.6)+
  geom_hline(yintercept = avg_score_diff, linetype="dashed")+
  geom_hline(yintercept = avg_score_diff+(avg_score_diff_SD), linetype="dashed", color="gray4")+
  geom_hline(yintercept = avg_score_diff-(avg_score_diff_SD), linetype="dashed", color="gray4")+
  xlab("Position")+
  ylab("Score Difference")+
  scale_x_continuous(breaks = seq(1059,1345, by=3))+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.background = element_blank())
plot(plot_met_ex_pos_difference)


plot_met_ex_pos_difference_2 <- ggplot()+
  geom_line(aes(x=avg_diff$pos, y=avg_diff$mean_diff), color = "#219ebc", size=0.6)+
  geom_hline(yintercept = 0)+ #no diff
  geom_hline(yintercept = avg_score_diff, linetype="dashed", color ="gray50")+
  geom_hline(yintercept = avg_score_diff+(avg_score_diff_SD), linetype="dashed", color="gray50")+
  geom_hline(yintercept = avg_score_diff-(avg_score_diff_SD), linetype="dashed", color="gray50")+
  xlab("Position")+
  ylab("Score Difference")+
  guides(x = "none")+
  scale_x_continuous(breaks = seq(1059,1345, by = 15))+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        axis.text = element_text(face="bold"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #axis.text.x=element_blank(),
        plot.background = element_blank())
plot(plot_met_ex_pos_difference_2)


JM_helix_diff = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N"&(pos %in% c(1062,1066, 1121,1125,1129)))
Chelix_diff = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1119,1122,1123,1126)))
ex14_Chelix_GOF = ex14_2SD_df_filtered %>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1119,1122,1123,1126)))
met_Chelix_GOF = met_2SD_df_filtered %>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1119,1122,1123,1126)))

plot_C_helix_corr<- ggplot()+
  geom_point(data =ex14_Chelix_GOF, aes(x=met_score, y=ex14_score), color = "#219ebc")+
  #geom_text(aes(x=Chelix_diff$met_score, y=Chelix_diff$ex14_score,label = Chelix_diff$hgvs))+
  scale_fill_manual(values=cbp1)+
  xlab("METdEx14") + ylab("MET")+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed",size=0.8)+
  geom_vline(xintercept = met_syn_mean, linetype="dashed",size=0.8)+
  ggtitle("C-helix")+
  ylim(-10,3)+
  xlim(-10,3)+
  theme_classic()+
  theme(text = element_text(size = 15, color="black"),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())
plot(plot_C_helix_corr)



plot_JM_helix_diff<- ggplot()+
  geom_boxplot(data=JM_helix_diff , aes(x=factor(pos), y= met_score),color="darkorange3", fill="orange", alpha=0.5)+
  geom_point(data=JM_helix_diff , aes(x=factor(pos), y= met_score),color="darkorange3" )+
  scale_fill_manual(values=cbp1)+
  ylim(-10,2)+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("METdEx14")+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())

plot_JM_helix_diff_2<- ggplot()+
  geom_boxplot(data =JM_helix_diff, aes(x=factor(pos), y=ex14_score), color="#219ebc", fill="#219ebc", alpha=0.3)+
  geom_point(data =JM_helix_diff, aes(x=factor(pos), y=ex14_score), color = "#219ebc")+
  scale_fill_manual(values=cbp1)+
  ylim(-10,2)+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("MET")+
  guides(y = "none")+
  theme_classic()+
  theme(text = element_text(size = 15, color="black"),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())
grid.arrange(plot_JM_helix_diff,plot_JM_helix_diff_2,ncol=2,widths=c(1.3,1))


plot_Chelix_diff<- ggplot()+
  geom_boxplot(data=Chelix_diff , aes(x=factor(pos), y= met_score),color="darkorange3", fill="orange", alpha=0.5)+
  geom_point(data=Chelix_diff , aes(x=factor(pos), y= met_score),color="darkorange3" )+
  scale_fill_manual(values=cbp1)+
  ylim(-10,2)+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("METdEx14")+
  geom_hline(yintercept = met_syn_mean, linetype="dashed",size=0.8)+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())

plot_Chelix_diff_2<- ggplot()+
  geom_boxplot(data =Chelix_diff, aes(x=factor(pos), y=ex14_score), color="#219ebc", fill="#219ebc", alpha=0.3)+
  geom_point(data =Chelix_diff, aes(x=factor(pos), y=ex14_score), color = "#219ebc")+
  scale_fill_manual(values=cbp1)+
  ylim(-10,2)+
  xlab("Position") + ylab("Activity Score")+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed",size=0.8)+
  ggtitle("MET")+
  guides(y = "none")+
  theme_classic()+
  theme(text = element_text(size = 15, color="black"),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())
grid.arrange(plot_Chelix_diff,plot_Chelix_diff_2,ncol=2,widths=c(1.3,1))

plot(plot_Chelix_diff_2)
#### helix-G correlations 

#G_helix_diff = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N"&(pos %in% c(1284,1287,1290,1292)))
G_helix_diff = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N"&(pos %in% c(1284,1287,1290,1291,1292,1293,
                                                                                                  1294,1297,1298,1299,1300)))

plot_G_helix_diff<- ggplot()+
  geom_point(data =G_helix_diff , aes(x=met_score, y=ex14_score), color = "#219ebc")+
  #geom_text(aes(x=G_helix_diff$met_score, y=G_helix_diff$ex14_score,label = G_helix_diff$hgvs))+
  scale_fill_manual(values=cbp1)+
  xlab("METdEx14") + ylab("MET")+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed",size=0.8)+
  geom_vline(xintercept = met_syn_mean, linetype="dashed",size=0.8)+
  ggtitle("G-helix loop")+
  ylim(-10,3)+
  xlim(-10,3)+
  theme_classic()+
  theme(text = element_text(size = 15, color="black"),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())
plot(plot_G_helix_diff)

plot_G_helix_diff_2<- ggplot()+
  geom_histogram(data =G_helix_diff , aes(x=met_score),fill = "#219ebc", alpha=0.4)+
  geom_histogram(data =G_helix_diff , aes(x=ex14_score), fill= "grey", alpha=0.4)+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed",size=0.8, color = "grey")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed",size=0.8, color = "#219ebc")+
  geom_vline(xintercept = 2*sd(met_syn$IL3_withdrawal_score), linetype="dashed",size=0.8, color = "blue")+
  geom_vline(xintercept = 2*sd(ex14_syn$IL3_withdrawal_score), linetype="dashed",size=0.8, color = "grey")+
  ggtitle("G-helix loop")+
  theme_classic()+
  theme(text = element_text(size = 15, color="black"),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())
plot(plot_G_helix_diff_2)

#### helix-D loop correlations 

D_helix_diff = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N"&(pos %in% c(1165, 1166, 1167, 1168, 1169, 1170,
                                                                                                  1171, 1172, 1173, 1174, 1175, 1176,
                                                                                                  1177, 1178, 1179)))

plot_D_helix_diff<- ggplot()+
  geom_point(data =D_helix_diff , aes(x=met_score, y=ex14_score), color = "mediumpurple3")+
  #geom_text(aes(x=D_helix_diff$met_score, y=D_helix_diff$ex14_score,label = D_helix_diff$hgvs))+
  scale_fill_manual(values=cbp1)+
  xlab("METdEx14") + ylab("MET")+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed",size=0.8)+
  geom_vline(xintercept = met_syn_mean, linetype="dashed",size=0.8)+
  ylim(-10,3)+
  xlim(-10,3)+
  theme_classic()+
  theme(text = element_text(size = 15, color="black"),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())
plot(plot_D_helix_diff)


#### JM/C helix correlations 

JM_C_helix_diff = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N"&(pos %in% c(1059, 1060, 1062, 1063, 1064, 1065, 1066, 1067, 1068,
                                                                                                     1069, 1070, 1117, 1118, 1119, 1120, 1121, 1122, 1123,
                                                                                                     1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132,
                                                                                                     1134)))

#JM_C_helix_diff = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N"&(pos %in% c(1068)))


plot_JM_C_helix_diff<- ggplot()+
  geom_point(data =JM_C_helix_diff, aes(x=met_score, y=ex14_score), color = "darkorange3")+
  geom_text(aes(x=JM_C_helix_diff$met_score, y=JM_C_helix_diff$ex14_score,label = JM_C_helix_diff$hgvs))+
  scale_fill_manual(values=cbp1)+
  xlab("METdEx14") + ylab("MET")+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed",size=0.8)+
  geom_vline(xintercept = met_syn_mean, linetype="dashed",size=0.8)+
  ggtitle("JM/C")+
  ylim(-10,3)+
  xlim(-10,3)+
  theme_classic()+
  theme(text = element_text(size = 15, color="black"),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank())
plot(plot_JM_C_helix_diff)

#####################################################################################################
#------------------------------- Comparisons of GOF and LOF characterization ------------------------
#####################################################################################################

# generate dataframe with GOF, LOF for met+/-ex14 scores +/-2 standard deviations from the mean of the respective synonymous distribution 
met_GOF_LOF_df <- subset(met_scores_ %>% filter(abs(IL3_withdrawal_score - mean(met_syn$IL3_withdrawal_score)) > 2*sd(met_syn$IL3_withdrawal_score)))
ex14_GOF_LOF_df <- subset(ex14_scores_ %>% filter(abs(IL3_withdrawal_score - mean(ex14_syn$IL3_withdrawal_score)) > 2*sd(ex14_syn$IL3_withdrawal_score)))

## what mutations are GOF and LOF for both MET+/-Ex14? ##
ex14_GOF_LOF<-data.frame(hgvs=ex14_GOF_LOF_df$hgvs,
                     ex14_SE = ex14_GOF_LOF_df$IL3_withdrawal_SE,
                     ex14_score=ex14_GOF_LOF_df$IL3_withdrawal_score,       
                     mutation_type=ex14_GOF_LOF_df$mutation_type,
                     pos=ex14_GOF_LOF_df$pos,
                     variants=ex14_GOF_LOF_df$variants)


met_GOF_LOF<-data.frame(hgvs=met_GOF_LOF_df$hgvs, 
                   met_SE = met_GOF_LOF_df$IL3_withdrawal_SE,
                   met_score =met_GOF_LOF_df$IL3_withdrawal_score,
                   mutation_type=met_GOF_LOF_df$mutation_type,
                   pos=met_GOF_LOF_df$pos, 
                   variants=met_GOF_LOF_df$variants)

GOF_LOF_compare <- merge(met_GOF_LOF, ex14_GOF_LOF)
GOF_LOF_compare$score_diff <- abs((GOF_LOF_compare$met_score)-(GOF_LOF_compare$ex14_score)) #difference 
GOF_LOF_compare$prop_error <- ((GOF_LOF_compare$met_SE)^2 +(GOF_LOF_compare$ex14_SE)^2)^(1/2) #propagated error
#GOF_LOF_compare <- GOF_LOF_compare %>% filter (GOF_LOF_compare$score_diff >= GOF_LOF_compare$prop_error )

#### GLOBAL GOF and LOF for MET+/-Ex14#### 

# plot the max score of GOF and LOF mutations on MET structures
met_GOF <- met_GOF_LOF_df %>% filter(IL3_withdrawal_score > 0) 
met_GOF_max <- met_GOF %>% group_by(pos) %>% summarise(met_GOF_max_scores= max(IL3_withdrawal_score, na.rm=TRUE))
met_GOF_max_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_GOF_max_3R7O, met_GOF_max, "met_GOF_max_scores")
write.pdb(x, file="met_GOF_max_3R7O.pdb")

met_LOF <- met_GOF_LOF_df %>% filter(IL3_withdrawal_score < 0) 
met_LOF_max <- met_LOF %>% group_by(pos) %>% summarise(met_LOF_max_scores= max(IL3_withdrawal_score, na.rm=TRUE))
met_LOF_max_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_LOF_max_3R7O, met_LOF_max, "met_LOF_max_scores")
write.pdb(x, file="met_LOF_max_3R7O.pdb")

met_GOF <- met_GOF_LOF_df %>% filter(IL3_withdrawal_score > 0) 
met_GOF_max <- met_GOF %>% group_by(pos) %>% summarise(met_GOF_max_scores= max(IL3_withdrawal_score, na.rm=TRUE))
met_GOF_max_2G15 <- read.pdb("2G15")
x = map_scores_pdb(met_GOF_max_2G15, met_GOF_max, "met_GOF_max_scores")
write.pdb(x, file="met_GOF_max_2G15.pdb")


# plot the max score of GOF and LOF mutations on +Ex14structures
ex14_GOF <- ex14_GOF_LOF_df %>% filter(IL3_withdrawal_score > 0) 
ex14_GOF_max <- ex14_GOF %>% group_by(pos) %>% summarise(ex14_GOF_max_scores= max(IL3_withdrawal_score, na.rm=TRUE))
ex14_GOF_max_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_GOF_max_3R7O, ex14_GOF_max, "ex14_GOF_max_scores")
write.pdb(x, file="ex14_GOF_max_3R7O.pdb")

met_GOF <- met_GOF_LOF_df %>% filter(IL3_withdrawal_score > 0) 
met_GOF_max <- met_GOF %>% group_by(pos) %>% summarise(met_GOF_max_scores= max(IL3_withdrawal_score, na.rm=TRUE))
met_GOF_max_2G15 <- read.pdb("2G15")
x = map_scores_pdb(met_GOF_max_2G15, met_GOF_max, "met_GOF_max_scores")
write.pdb(x, file="met_GOF_max_2G15.pdb")

ex14_LOF <- ex14_GOF_LOF_df %>% filter(IL3_withdrawal_score < -0.1426368) 
ex14_LOF_avg <- ex14_LOF %>% group_by(pos) %>% summarise(ex14_LOF_avg_scores= mean(IL3_withdrawal_score, na.rm=TRUE))
ex14_LOF_avg_2G15 <- read.pdb("2G15")
x = map_scores_pdb(ex14_LOF_avg_2G15, ex14_LOF_avg, "ex14_LOF_avg_scores")
write.pdb(x, file="ex14_LOF_avg_2G15.pdb")
ex14_LOF <- ex14_GOF_LOF_df %>% filter(IL3_withdrawal_score < -0.1426368) 

ex14_LOF_avg <- ex14_LOF %>% group_by(pos) %>% summarise(ex14_LOF_avg_scores= abs(sum(IL3_withdrawal_score, na.rm=TRUE)))
ex14_LOF_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_LOF_avg_3R7O, ex14_LOF_avg, "ex14_LOF_avg_scores")
write.pdb(x, file="ex14_LOF_avg_3R7O.pdb")

ex14_LOF <- ex14_GOF_LOF_df %>% filter(IL3_withdrawal_score < -0.1426368) 
ex14_LOF_avg <- ex14_LOF %>% group_by(pos) %>% summarise(ex14_LOF_avg_scores= abs(IL3_withdrawal_score, na.rm=TRUE))
ex14_LOF_avg_3RHK <- read.pdb("3RHK")
x = map_scores_pdb(ex14_LOF_avg_3RHK, ex14_LOF_avg, "ex14_LOF_avg_scores")
write.pdb(x, file="ex14_LOF_avg_3RHK.pdb")

##### Direct positional comparisons #######

# plot positions that are LOF for both constructs 
LOF_compare <- subset(GOF_LOF_compare, met_score < 0 & ex14_score < 0)
GOF_compare <- subset(GOF_LOF_compare, met_score > 0 & ex14_score > 0 )

LOF_compare_met_abs_sum <- LOF_compare %>% group_by(pos) %>% summarise(met_LOF_compare = abs(sum(met_score, na.rm=TRUE)))
LOF_compare_ex14_abs_sum <- LOF_compare %>% group_by(pos) %>% summarise(ex14_LOF_compare = abs(sum(ex14_score, na.rm=TRUE)))

met_LOF_compare_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_LOF_compare_3R7O, LOF_compare_met_abs_sum, "met_LOF_compare")
write.pdb(x, file="met_LOF_compare_3R7O.pdb")
ex14_LOF_compare_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_LOF_compare_3R7O, LOF_compare_ex14_abs_sum, "ex14_LOF_compare")
write.pdb(x, file="ex14_LOF_compare_3R7O.pdb")

LOF_compare_met_abs_sum <- LOF_compare %>% group_by(pos) %>% summarise(met_LOF_compare = abs(sum(met_score, na.rm=TRUE)))
LOF_compare_ex14_abs_sum <- LOF_compare %>% group_by(pos) %>% summarise(ex14_LOF_compare = abs(sum(ex14_score, na.rm=TRUE)))

met_GOF_compare_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_GOF_compare_3R7O, GOF_compare_met_abs_sum, "met_GOF_compare")
write.pdb(x, file="met_GOF_compare_3R7O.pdb")

ex14_LOF_compare_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_GOF_compare_3R7O, LOF_compare_ex14_abs_sum, "ex14_GOF_compare")
write.pdb(x, file="ex14_GOF_compare_3R7O.pdb")

met_GOF_compare_2G15 <- read.pdb("2G15")
x = map_scores_pdb(met_GOF_compare_2G15, GOF_compare_met_abs_sum, "met_GOF_compare")
write.pdb(x, file="met_GOF_compare_2G15.pdb")
ex14_GOF_compare_2G15 <- read.pdb("2G15")
x = map_scores_pdb(ex14_GOF_compare_2G15, GOF_compare_ex14_abs_sum, "ex14_GOF_compare")
write.pdb(x, file="ex14_GOF_compare_2G15.pdb")


#####################################################################################################
### -------------------------------- center of mass cluster of GOF  ----------------------------
#####################################################################################################

# inspired by the SRC validation of GOF in 3D space 
# utilizes the center of mass of each amino acid from the bio3d package 
# hierarchical clustering with dendextend from the R packages 
# using GOF mutations -- GOF = +2.5 SD from mean of synonymous distirbution 

# calculates the absolute sum of the score for each position 
met_GOF_pos <-  met_GOF_LOF_df %>% filter(IL3_withdrawal_score > 0) %>%
  group_by(pos) %>% summarise(met_GOF_ = abs(sum(IL3_withdrawal_score, na.rm=TRUE)))


### gets centroid for each residue and appends it to the corresponding position and abs(sum(score))
pdb_active <- read.pdb("3R7O") # open PDB
met_GOF_residue_pos <-c(1059,1073,1075,1078,1096,1098,1112,1144,1146,
                        1148,1157,1159,1160,1161,1163,1170,1177,1179,1180,1190,
                        1198,1217,1227,1228,1233,1249,1261,1278,1298,1310,1314,
                        1324,1326,1334,1335,1336)


met_GOF_centroids <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("pos", "x", "y","z"))))
for (i in met_GOF_residue_pos){
  inds <- atom.select(pdb_active, chain="A",resno=i)
  met_GOF_centroids  <- rbind(met_GOF_centroids,c(i,com(pdb_active,inds)))
  print(com(pdb_active, inds))
  print(i)
}

colnames(met_GOF_centroids ) <- c("pos","x","y","z")

### plots dendrogram
d <- dist(met_GOF_centroids[,-1])
hc <- hclust(d)
plot(hc,lab=met_GOF_centroids$pos,hang = -1, cex = 0.6)

hc$labels = met_GOF_centroids$pos
dend <- as.dendrogram(hc)
d1=color_branches(dend,k=5)
plot(d1) # selective coloring of branches :)
d2=color_branches(d1,k=5) # auto-coloring 5 clusters of branches.
plot(d2)

######### Unbiased CLustering of GOF mutaitons ######

ex14_GOF_pos <-  ex14_GOF_LOF_df %>% filter(IL3_withdrawal_score > 0) %>%
  group_by(pos) %>% summarise(ex14_GOF_ = abs(sum(IL3_withdrawal_score, na.rm=TRUE)))
ex14_GOF_residue_pos <-c(1060,1077,1078,1096,1106,1122,1123,1126,1135,1136,1138,1146,
                         1150,1159,1161,1165,1166,1170,1174,1177,1178,1180,1182,1184,
                         1188,1190,1197,1207,1211,1215,1217,1233,1236,1237,1239,1244,
                         1293,1294,1299,1303,1304,1311,1321,1324,1326,1328,
                         1334,1341)


ex14_GOF_centroids <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("pos", "x", "y","z"))))
for (i in ex14_GOF_residue_pos){
  inds2 <- atom.select(pdb_active, chain="A",resno=i)
  ex14_GOF_centroids  <- rbind(ex14_GOF_centroids,c(i,com(pdb_active,inds2)))
  print(com(pdb_active, inds2))
  print(i)
}

colnames(ex14_GOF_centroids ) <- c("pos","x","y","z")

### plots dendrogram
dd <- dist(ex14_GOF_centroids[,-1])
hh <- hclust(dd)
plot(hh,lab=ex14_GOF_centroids$pos,hang = -1, cex = 0.6)

hh$labels = ex14_GOF_centroids$pos
dend <- as.dendrogram(hh)
d1=color_branches(dend,k=5)
plot(d1) # selective coloring of branches :)
d2=color_branches(d1,k=5) # auto-coloring 5 clusters of branches.
plot(d2)


######### Unbiased Clustering of LOF  mutaitons ######

met_LOF_pos <-  met_GOF_LOF_df %>% filter(IL3_withdrawal_score < 0) %>%
  group_by(pos) %>% summarise(met_LOF_ = mean(IL3_withdrawal_score, na.rm=TRUE))

pdb_active <- read.pdb("3R7O") # open PDB
### these are ther positions w/ the greatest magnitude of effect at a position -- top 40 positions 
#met_LOF_residue_pos <-c(1265,1267,1110,1191,1128,1287,1275,1220,1301,1280,1227,
                        #1204,1246,1253,1235,1206,1225,1221,1274,1108,1090,1254,
                        #1087,1095,1208,1327,1203,1234,1209,1270,1202,1153,1251,
                        #1229,1222,1284,1165,1283,1195,1271)
met_LOF_residue_pos <-c(1282,1291,1257,1287,1328,1313,1208,1089,1304,1191,1097,1206,
                        1265,1296,1220,1280,1305,1108,1300,1226,1164,1254,1253,1095,
                        1182,1292,1295,1252,1110,1290,1218,1308,1245,1316,1125,1156,
                        1128,1225,1234,1267)


met_LOF_centroids <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("pos", "x", "y","z"))))
for (i in met_LOF_residue_pos){
  inds <- atom.select(pdb_active, chain="A",resno=i)
  met_LOF_centroids  <- rbind(met_LOF_centroids,c(i,com(pdb_active,inds)))
  print(com(pdb_active, inds))
  print(i)
}

colnames(met_LOF_centroids ) <- c("pos","x","y","z")

### plots dendrogram
d <- dist(met_LOF_centroids[,-1])
hc <- hclust(d)
plot(hc,lab=met_LOF_centroids$pos,hang = -1, cex = 0.6)

hc$labels = met_LOF_centroids$pos
dend <- as.dendrogram(hc)
d1=color_branches(dend,k=5)
plot(d1) # selective coloring of branches :)
d2=color_branches(d1,k=5) # auto-coloring 5 clusters of branches.
plot(d2)



ex14_LOF_pos <-  ex14_GOF_LOF_df %>% filter(IL3_withdrawal_score < 0) %>%
  group_by(pos) %>% summarise(ex14_LOF_ = mean(IL3_withdrawal_score, na.rm=TRUE))

pdb_active <- read.pdb("3R7O") # open PDB
### these are ther positions w/ the greatest magnitude of effect at a position -- top 40 positions abs(sum(score))
#ex14_LOF_residue_pos <-c(1301,1265,1284,1253,1235,1267,1287,1191,1246,1227,1204,
                         #1206,1300,1225,1087,1153,1110,1127,1071,1222,1234,1203,
                         #1275,1260,1221,1290,1128,1248,1092,1090,1208,1296,1209,
                         #1274,1223,1254,1125,1108,1124,1164)

### these are ther positions w/ the greatest magnitude of effect at a position -- top 40 positions mean
ex14_LOF_residue_pos <-c(1180,1064,1217,1136,1063,1071,1208,1253,1301,1065,1265,
                         1229,1165,1284,1256,1092,1271,1298,1235,1186,1267,1156,
                         1312,1213,1287,1264,1280,1123,1169,1206,1300,1283,1090,
                         1225,1296,1120,1260,1246,1129,1249)



ex14_LOF_centroids <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("pos", "x", "y","z"))))
for (i in ex14_LOF_residue_pos){
  inds <- atom.select(pdb_active, chain="A",resno=i)
  ex14_LOF_centroids  <- rbind(ex14_LOF_centroids,c(i,com(pdb_active,inds)))
  print(com(pdb_active, inds))
  print(i)
}

colnames(ex14_LOF_centroids ) <- c("pos","x","y","z")

### plots dendrogram
d <- dist(ex14_LOF_centroids[,-1])
hc <- hclust(d)
plot(hc,lab=ex14_LOF_centroids$pos,hang = -1, cex = 0.6)

hc$labels = ex14_LOF_centroids$pos
dend <- as.dendrogram(hc)
d1=color_branches(dend,k=5)
plot(d1) # selective coloring of branches :)
d2=color_branches(d1,k=5) # auto-coloring 5 clusters of branches.
plot(d2)

#####################################################################################################
###----------------------------positional scores scatter plot --------------------------------------
#####################################################################################################

plot <- ggplot(data=merged_datasets) + 
  geom_point(shape=21,aes(x=pos,y=met_score))+
  xlab("MET Position") + ylab("MET fitness score")+
  ylim(-15,15)
plot(plot)

plot <- ggplot(data=merged_datasets) + 
  geom_point(shape=21,aes(x=pos,y=ex14_score))+
  xlab("Ex14 Position") + ylab("Ex24 fitness score")+
  ylim(-15,15)
plot(plot)

#####################################################################################################
###------- plot SE propagation of error as a histogram
#####################################################################################################
plot<-ggplot(data=merged_datasets, aes(x=prop_error)) + 
  geom_histogram(aes(y =..density..), 
                 colour = "black", 
                 fill = "#56B4E9",
                 alpha=.2) +
  xlim(-1,4)+
  geom_vline(xintercept=c(1.2), linetype="dotted", color ='red')
  #stat_function(fun = dnorm, args = list(mean = mean(merged_datasets$prop_error), sd = sd(merged_datasets$prop_error)))
plot(plot)


#####################################################################################################
### averages as scatter plot with error bars 
#####################################################################################################
plot <- ggplot(data = merged_datasets_avg, aes(x=met_avg_IL3_withdrawal, y=ex14_avg_IL3_withdrawal)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ex14_avg_IL3_withdrawal-ex14_avg_IL3_withdrawal_SE, 
                    ymax = ex14_avg_IL3_withdrawal+ex14_avg_IL3_withdrawal_SE)) + 
  geom_errorbarh(aes(xmin = met_avg_IL3_withdrawal-met_avg_IL3_withdrawal_SE, 
                     xmax = met_avg_IL3_withdrawal+met_avg_IL3_withdrawal_SE))
plot(plot)

#####################################################################################################
###---------------------------- Plot average fitness per POSTITION  ---------------------------------
#####################################################################################################

plot <- ggplot(merged_datasets_avg, aes(x=met_avg_IL3_withdrawal, y=ex14_avg_IL3_withdrawal))+
  geom_point(shape=21)+
  geom_smooth(method=lm, color='black', fill='orange')
  theme_linedraw()+
  xlab("MET-Ex14") + ylab("MET+Ex14")
plot(plot)

plot <- ggplot(merged_datasets_avg, aes(x=met_avg_IL3, y=ex14_avg_IL3))+
  geom_point(shape=21)+
  #geom_density_2d()+
  theme_linedraw()+
  xlab("MET-Ex14") + ylab("MET+Ex14")
plot(plot)

plot <- ggplot(data = merged_datasets_avg, aes(x=met_avg_IL3_withdrawal, y=ex14_avg_IL3_withdrawal, label = pos)) +
  geom_point() + 
  #geom_errorbar(aes(ymin = ex14_avg_IL3_withdrawal-ex14_avg_IL3_withdrawal_SE, 
                    #ymax = ex14_avg_IL3_withdrawal+ex14_avg_IL3_withdrawal_SE)) + 
  #geom_errorbarh(aes(xmin = met_avg_IL3_withdrawal-met_avg_IL3_withdrawal_SE, 
                     #xmax = met_avg_IL3_withdrawal+met_avg_IL3_withdrawal_SE))+
  geom_text()
plot(plot)



#####################################################################################################
###----------- Quadrant analysis of statitically different mutations ------------------------------
#####################################################################################################

# IL3 withdrawal
#df_quad1<-subset(merged_datasets, ex14_score > -5 & met_score > -5) # quadarant 1
#df_quad2<-subset(merged_datasets, ex14_score > -5 & met_score < -5) # quadarant 2
#df_quad3<-subset(merged_datasets, ex14_score < -5 & met_score < -5) # quadarant 3
#df_quad4<-subset(merged_datasets, ex14_score < -5 & met_score > -5 & met_score < -1) # quadarant 4


df_quad1<-subset(scores_filtered_diff, ex14_score > -5 & met_score > -5) # quadarant 1
df_quad2<-subset(scores_filtered_diff, ex14_score > -5 & met_score < -5) # quadarant 2
df_quad3<-subset(scores_filtered_diff, ex14_score < -5 & met_score < -5) # quadarant 3
df_quad4<-subset(scores_filtered_diff, ex14_score < -2.5 & met_score > -2.5) # quadarant 4
df_quad5<-subset(scores_filtered_diff, ex14_score > 0 & met_score > 0)




quad1_avg <- df_quad1 %>% group_by(pos) %>% summarise(quad1_met = mean(met_score, na.rm=TRUE),
                                                      quad1_ex14 = mean(ex14_score, na.rm=TRUE))

quad2_avg <- df_quad2 %>% group_by(pos) %>% summarise(quad2_met = mean(met_score, na.rm=TRUE),
                                                      quad2_ex14 = mean(ex14_score, na.rm=TRUE))

quad3_avg <- df_quad3 %>% group_by(pos) %>% summarise(quad3_met = mean(met_score, na.rm=TRUE),
                                                      quad3_ex14 = mean(ex14_score, na.rm=TRUE))

quad4_avg <- df_quad4 %>% group_by(pos) %>% summarise(quad4_met = mean(met_score, na.rm=TRUE),
                                                      quad4_ex14 = mean(ex14_score, na.rm=TRUE))

quad4_max <- df_quad4 %>% group_by(pos) %>% summarise(quad4_met = max(met_score, na.rm=TRUE),
                                                      quad4_ex14 = max(ex14_score, na.rm=TRUE))

quad4_abs_sum <- df_quad4 %>% group_by(pos) %>% summarise(quad4_met = sum(met_abs_score, na.rm=TRUE),
                                                      quad4_ex14 = sum(ex14_abs_score, na.rm=TRUE))
quad4_abs_sum$diff <-abs((quad4_abs_sum$quad4_met)-(quad4_abs_sum$quad4_ex14)) #absolute difference 



quad5_avg <- df_quad5 %>% group_by(pos) %>% summarise(quad5_met = mean(met_score, na.rm=TRUE),
                                                      quad5_ex14 = mean(ex14_score, na.rm=TRUE))


###--------------------- quad4 mutations by phiochem -----------------------------
df_quad4_extreme<-subset(merged_datasets, ex14_score < -5 & met_score > 2)
quad4_hydrophobic = df_quad4_extreme %>% filter(variants %in% c("V","L","I","M","P","G","A","W","F"))
plot <- ggplot(data=quad4_hydrophobic, aes(x=met_score, y=ex14_score, color = mutation_type, label=hgvs)) + 
  geom_point(shape = 1, size = 3) +
  theme_pubr() +
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("Hydrophobic")+
  geom_text( color = "black",position = position_nudge(y = 0.2))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
plot(plot)


df_quad5_extreme<-subset(merged_datasets, ex14_score < -3 & met_score > 2)
quad5_aromatic = df_quad5_extreme %>% filter(variants %in% c("Y","W","F"))
plot <- ggplot(data=quad5_aromatic, aes(x=met_score, y=ex14_score, color = mutation_type, label=hgvs)) + 
  geom_point(shape = 1, size = 3) +
  theme_pubr() +
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("Aromatic")+
  geom_text( color = "black",position = position_nudge(y = 0.2))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
plot(plot)


df_quad1_extreme<-subset(merged_datasets, ex14_score >1 & met_score > -5)
plot <- ggplot(data=quad1_acidic, aes(x=met_score, y=ex14_score, color = mutation_type, label=hgvs)) + 
  geom_point(shape = 1, size = 3) +
  theme_pubr() +
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("acidic")+
  geom_text( color = "black",position = position_nudge(y = 0.2))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
plot(plot)

#####################################################################################################
###----------- map avg positional mutation scores from each quadrant on a pdb  ----------------------
#####################################################################################################

### color with: spectrum b, firebrick_white_sky_blue, minimum = -5, maximum=5

### quad 1 pdb
ex14_quad1_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_quad1_avg_3R7O, quad1_avg , "quad1_ex14")
write.pdb(x, file="ex14_quad1_avg_3R7O.pdb")

met_quad1_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_quad1_avg_3R7O, quad1_avg, "quad1_met")
write.pdb(x, file="met_quad1_avg_3R7O..pdb")


### quad 2 pdb
ex14_quad2_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_quad2_avg_3R7O, quad2_avg , "quad2_ex14")
write.pdb(x, file="ex14_quad2_avg_3R7O.pdb")

met_quad2_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_quad2_avg_3R7O, quad2_avg, "quad2_met")
write.pdb(x, file="met_quad2_avg_3R7O.pdb")


### quad 3 pdb
ex14_quad3_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_quad3_avg_3R7O, quad3_avg , "quad3_ex14")
write.pdb(x, file="ex14_quad3_avg_3R7O.pdb")

met_quad3_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_quad3_avg_3R7O, quad3_avg, "quad3_met")
write.pdb(x, file="met_quad3_avg_3R7O.pdb")


### quad 4 pdb
ex14_quad4_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_quad4_avg_3R7O, quad4_avg , "quad4_ex14")
write.pdb(x, file="ex14_quad4_avg_3R7O.pdb")

met_quad4_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_quad4_avg_3R7O, quad4_avg, "quad4_met")
write.pdb(x, file="met_quad4_avg_3R7O.pdb")

### quad 4 pdb
ex14_quad4_max_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_quad4_max_3R7O, quad4_max , "quad4_ex14")
write.pdb(x, file="ex14_quad4_max_3R7O.pdb")

met_quad4_max_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_quad4_max_3R7O, quad4_max, "quad4_met")
write.pdb(x, file="met_quad4_max_3R7O.pdb")

### quad 4 pdb
ex14_quad4_abs_sum_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_quad4_abs_sum_3R7O, quad4_abs_sum , "quad4_ex14")
write.pdb(x, file="ex14_quad4_abs_sum_3R7O.pdb")

met_quad4_abs_sum_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_quad4_abs_sum_3R7O, quad4_abs_sum, "quad4_met")
write.pdb(x, file="met_quad4_abs_sum_3R7O.pdb")


quad4_abs_sum_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(quad4_abs_sum_3R7O, quad4_abs_sum, "diff")
write.pdb(x, file="quad4_abs_sum_3R7O.pdb")


### quad 5 pdb
ex14_quad5_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_quad5_avg_3R7O, quad5_avg , "quad5_ex14")
write.pdb(x, file="ex14_quad5_avg_3R7O.pdb")

met_quad5_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_quad5_avg_3R7O, quad5_avg, "quad5_met")
write.pdb(x, file="met_quad5_avg_3R7O..pdb")

ex14_3R7O<- read.pdb("2G15")
x = map_scores_pdb(scores_filtered_diff, ex14_scores, "avg_IL3")
write.pdb(x, file="ex14_2G15_avg_IL3.pdb")


#####################################################################################################
##----------------------------- Cancer and Resistance Mutations  ------------------------------------
#####################################################################################################

###---------------- plot known MET KD cancer mutations as distributions ----------------------------


# cancer positions

met_cancer_mut = met_scores_  %>% filter((hgvs %in% c("p.(M192T)","p.(V12E)","p.(Y172H)","p.(V12M)","p.(V12G)","p.(F142I)",
                                                      "p.(H36Y)","p.(Y172C)","p.(D170Y)","p.(D170A)","p.(Y177H)","p.(N42D)",
                                                      "p.(V34I)","p.(D170H)","p.(Q65R)","p.(H104R)","p.(D228G)","p.(L254V)",
                                                      "p.(M134I)","p.(T38S)","p.(F66L)","p.(R90Q)","p.(H180L)","p.(E267K)",
                                                      "p.(L147V)","p.(D122Y)","p.(L128F)","p.(V229I)","p.(R269C)","p.(A193T)",
                                                      "p.(V5A)","p.(P15S)","p.(L128R)","p.(N80K)","p.(L128I)","p.(K141T)",
                                                      "p.(E256G)","p.(R221I)","p.(C250S)","p.(C33Y)","p.(L82F)","p.(P253L)",
                                                      "p.(R221K)","p.(M102I)","p.(H116N)","p.(P95Q)","p.(Q200K)","p.(P2S)",
                                                      "p.(D41G)","p.(V63A)","p.(S178R)","p.(G105R)")))

ex14_cancer_mut = ex14_scores_ %>% filter((hgvs %in% c("p.(M192T)","p.(V12E)","p.(Y172H)","p.(V12M)","p.(V12G)","p.(F142I)",
                                                       "p.(H36Y)","p.(Y172C)","p.(D170Y)","p.(D170A)","p.(Y177H)","p.(N42D)",
                                                       "p.(V34I)","p.(D170H)","p.(Q65R)","p.(H104R)","p.(D228G)","p.(L254V)",
                                                       "p.(M134I)","p.(T38S)","p.(F66L)","p.(R90Q)","p.(H180L)","p.(E267K)",
                                                       "p.(L147V)","p.(D122Y)","p.(L128F)","p.(V229I)","p.(R269C)","p.(A193T)",
                                                       "p.(V5A)","p.(P15S)","p.(L128R)","p.(N80K)","p.(L128I)","p.(K141T)",
                                                       "p.(E256G)","p.(R221I)","p.(C250S)","p.(C33Y)","p.(L82F)","p.(P253L)",
                                                       "p.(R221K)","p.(M102I)","p.(H116N)","p.(P95Q)","p.(Q200K)","p.(P2S)",
                                                       "p.(D41G)","p.(V63A)","p.(S178R)","p.(G105R)")))

met_cancer_df<-data.frame(hgvs=met_cancer_mut$hgvs,
                   met_score=met_cancer_mut$IL3_withdrawal_score,
                   mutation_type=met_cancer_mut$mutation_type,
                   pos=met_cancer_mut$pos)

ex14_cancer_df<-data.frame(hgvs=ex14_cancer_mut$hgvs,
                     ex14_score=ex14_cancer_mut$IL3_withdrawal_score,       
                     mutation_type=ex14_cancer_mut$mutation_type, 
                     pos=ex14_cancer_mut$pos)

merged_cancer_datasets<-merge(met_cancer_df,ex14_cancer_df)


plot_known_met_cancer <- ggplot(met_cancer_mut) +
  geom_histogram(aes(x = IL3_withdrawal_score),color="black",fill="deepskyblue2", alpha=0.5) +
  xlim(-10,10)+
  xlab("Activity score, METdEx14 KD cancer-associated mutations") + 
  ylab("Count") + 
  ggtitle("METdEx14 KD cancer-associated mutations")+
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_classic()
plot(plot_known_met_cancer)


plot <- ggplot(data=merged_cancer_datasets, aes(x=met_score)) + 
  geom_point(shape=21,aes(y=ex14_score))+
  theme_linedraw()+
  xlab("MET-Ex14") + ylab("MET+Ex14")
plot(plot)

plot_mcanmut <- ggplot(data=merged_cancer_datasets, aes(x=pos,y=met_score, label=hgvs)) + 
  geom_point(shape=21)+
  ylab("MetdEx14 cancer mutation")+
  geom_text(aes(label=hgvs),hjust=0, vjust=0,size=3)+
  theme_linedraw()
plot(plot_mcanmut)

plot_excanmut <- ggplot(data=merged_cancer_datasets, aes(x=pos,y=ex14_score, label=hgvs)) + 
  geom_point(shape=21)+
  ylab("Met+Ex14 cancer mutation")+
  geom_text(aes(label=hgvs),hjust=0, vjust=0,size=3)+
  theme_linedraw()
plot(plot_excanmut)


plot <- ggplot(merged_datasets, aes(x=met_score, y=ex14_score))+
  geom_point(shape=21,aes(colour=mutation_type))+
  #geom_density_2d()+
  theme_linedraw()+
  xlab("MET-Ex14") + ylab("MET+Ex14")
plot(plot)


ggsave("met_cancer_mut_DFE.pdf", height = 4, width = 5)
ggsave("met_cancer_mut_DFE.png", height = 4, width = 5)

plot <- ggplot(ex14_cancer_mut) +
  geom_histogram(aes(x = IL3_withdrawal_score),color="black",fill="deepskyblue2", alpha=0.5) +
  xlim(-10,10)+
  xlab("Activity score, Ex14 KD cancer-associated mutations") + 
  ylab("Count") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot)

ggsave("ex14_cancer_mut_DFE.pdf", height = 4, width = 5)
ggsave("ex14_cancer_mut_DFE.png", height = 4, width = 5)


#####################################################################################################
###---------------- plot MET cancer mutation POSITIONS as distributions -----------------------------
#####################################################################################################

### all associated cancer positions from cBioPortal 
met_cancer_pos = met_scores_  %>% filter((pos %in% c(1060,1063,1070,1073,1091,1092,1094,1096,1099,
                                                    1100,1121,1123,1124,1138,1140,1148,1153,1160,1163,
                                                    1174,1180,1186,1192,1199,1200,1205,1228,1230,1235,1236,
                                                    1238,1250,1251,1258,1279,1286,1287,1308,1311,1312,1314,
                                                    1319,1325,1327)))

ex14_cancer_pos = ex14_scores_ %>% filter((pos %in% c(1060,1063,1070,1073,1091,1092,1094,1096,1099,
                                                     1100,1121,1123,1124,1138,1140,1148,1153,1160,1163,
                                                     1174,1180,1186,1192,1199,1200,1205,1228,1230,1235,1236,
                                                     1238,1250,1251,1258,1279,1286,1287,1308,1311,1312,1314,
                                                     1319,1325,1327)))

met_cancer_pos_df<-data.frame(hgvs=met_cancer_pos$hgvs,
                          met_score=met_cancer_pos$IL3_withdrawal_score,
                          mutation_type=met_cancer_pos$mutation_type,
                          pos=met_cancer_pos$pos)

ex14_cancer_pos_df<-data.frame(hgvs=ex14_cancer_pos$hgvs,
                           ex14_score=ex14_cancer_pos$IL3_withdrawal_score,       
                           mutation_type=ex14_cancer_pos$mutation_type, 
                           pos=ex14_cancer_pos$pos)

merged_cancer_pos_datasets<-merge(met_cancer_pos_df,ex14_cancer_pos_df)

## plots raw cancer-associated position scores
# met_syn_mean = -0.4341432
# ex14_syn_mean=-0.1426368

plot_mcan <- ggplot(met_cancer_pos) +
  geom_histogram(aes(x = IL3_withdrawal_score),color="black",fill="darkseagreen3", alpha=0.5) +
  xlab("Acivity Score") + 
  ylab("Count") + 
  xlim(-10,3)+
  geom_vline(xintercept=met_syn_mean, linetype="dotted",color="red")+
  geom_text(aes(x=met_syn_mean-4, label="mean = -0.43", y=125), colour="red", angle=0)+
  ggtitle("METdEx14 cancer positions")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())


plot_ecan <- ggplot(ex14_cancer_pos) +
  geom_histogram(aes(x = IL3_withdrawal_score),color="black",fill="darkseagreen3", alpha=0.5) +
  xlab("Activity score") + 
  xlim(-10,3)+
  geom_vline(xintercept=ex14_syn_mean, linetype="dotted",color="red")+
  geom_text(aes(x=ex14_syn_mean-4, label="mean = -0.14", y=125), colour="red", angle=0)+
  ggtitle("MET+Ex14 cancer positions")+
  ylab("Count") +
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

#### plot cancer mutations that are 2 SD above WT mean and have low error 

met_cancer_2SD <- subset(met_cancer_pos %>% 
                           filter(
                             (abs(IL3_withdrawal_score - mean(met_syn$IL3_withdrawal_score)) > 2*sd(met_syn$IL3_withdrawal_score)) 
                             & (IL3_withdrawal_SE <= met_wt_SE_mean))
                         )

plot_mcan_scatter <- ggplot(met_cancer_2SD %>% filter(IL3_withdrawal_score>0),aes(x=pos,y=IL3_withdrawal_score, label=hgvs)) +
  geom_point()+
  ggtitle("METdEx14 GOF cancer pos")+
  geom_text(hjust=0, vjust=0)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
  

ex14_cancer_2SD <- subset(ex14_cancer_pos %>% 
                           filter(
                             (abs(IL3_withdrawal_score - mean(ex14_syn$IL3_withdrawal_score)) > 2*sd(ex14_syn$IL3_withdrawal_score)) 
                             & (IL3_withdrawal_SE <= ex14_wt_SE_mean))
                          )

plot_ecan_scatter <- ggplot(ex14_cancer_2SD%>% filter(IL3_withdrawal_score>0),aes(x=pos,y=IL3_withdrawal_score, label=hgvs)) +
  geom_point()+
  ggtitle("MET+Ex14 GOF cancer pos")+
  geom_text(hjust=0, vjust=0)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())


plot_grid(plot_ecan, plot_ecan_scatter,plot_mcan,plot_mcan_scatter, ncol=2,nrow=2)

### plots distributions of cancer mutations 

cancer_color_pallet <- c("#999999","#03C03C")

data1 = met_cancer_mut$IL3_withdrawal_score
data2 = met_missense$IL3_withdrawal_score
d = data.frame(x = c(data1,data2), 
               type=rep(c("Known cancer mutations", "All Mutations"), 
                        c(length(data1), length(data2))))

plot_met_cancer_1 <- ggplot(d, aes(x, color=type, fill=type)) +
  scale_colour_manual(values=cancer_color_pallet)+
  scale_fill_manual(values=cancer_color_pallet)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.5,
                 mapping = aes(y = stat(ncount)))+
  geom_vline(xintercept = ex14_syn_mean,linetype="dashed")+
  ggtitle("METdEx14 cancer")+
  xlab("Activity Score")+
  ylab("Normalizaed Count")+
  xlim(-10,3)+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")



data3 = ex14_cancer_mut$IL3_withdrawal_score
data4 = ex14_missense$IL3_withdrawal_score
c = data.frame(x = c(data3,data4), 
               type=rep(c("Known cancer mutations", "All mutations"), 
                        c(length(data3), length(data4))))

plot_ex14_cancer_1 <- ggplot(c, aes(x, color=type, fill=type)) +
  scale_colour_manual(values=cancer_color_pallet)+
  scale_fill_manual(values=cancer_color_pallet)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.5,
                 mapping = aes(y = stat(ncount)))+
  geom_vline(xintercept = met_syn_mean,linetype="dashed")+
  ggtitle("MET+Ex14 cancer ")+
  xlab("Activity Score")+
  ylab("Normalizaed Count")+
  xlim(-10,3)+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1.7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")

plot_grid(plot_ex14_cancer_1,plot_met_cancer_1, ncol=2,nrow=1)


#####################################################################################################
###---------------- plot Cancer Frequencies Detected in Patients  ----------------------------------
#####################################################################################################

cancer_freq <-data.frame(
  hgvs = c("p.(M192T)","p.(V12E)","p.(Y172H)","p.(V12M)","p.(V12G)","p.(F142I)",
           "p.(H36Y)","p.(Y172C)","p.(D170Y)","p.(D170A)","p.(Y177H)","p.(N42D)",
           "p.(V34I)","p.(D170H)","p.(Q65R)","p.(H104R)","p.(D228G)","p.(L254V)",
           "p.(M134I)","p.(T38S)","p.(F66L)","p.(R90Q)","p.(H180L)","p.(E267K)",
           "p.(L147V)","p.(D122Y)","p.(L128F)","p.(V229I)","p.(R269C)","p.(A193T)",
           "p.(V5A)","p.(P15S)","p.(L128R)","p.(N80K)","p.(L128I)","p.(K141T)",
           "p.(E256G)","p.(R221I)","p.(C250S)","p.(C33Y)","p.(L82F)","p.(P253L)",
           "p.(R221K)","p.(M102I)","p.(H116N)","p.(P95Q)","p.(Q200K)","p.(P2S)",
           "p.(D41G)","p.(V63A)","p.(S178R)"),
  protein_change = c("M1250T","V1070E","Y1230H","V1070M","V1070G","F1200I",
                     "H1094Y","Y1230C","D1228Y","D1228A","Y1235H","N1100D",
                     "V1092I","D1228H","Q1123R","H1162R","D1286G","L1312V",
                     "M1192I","T1096S","F1124L","R1148Q","H1238L","E1325K",
                     "L1205V","D1180Y","L1186F","V1287I","R1327C","A1251T",
                     "V1063A","P1073S","L1186R","N1138K","L1186I","K1199T",
                     "E1314G","R1279I","C1308S","C1091Y","L1140F","P1311L",
                     "R1279K","M1160I","H1174N","P1153Q","Q1258K","P1060S",
                     "D1099G","V1121A","S1236R"),
  cBioPortal_frq = as.numeric(c(118,182,206,112,112,37,243,79,15,5,11,33,155,24,
                                7,98,44,11,19,17,9,167,53,65,90,14,32,13,3,61,34,
                                16,10,17,11,14,21,17,23,4,17,72,115,59,4,4,4,114,
                                51,236,55))
  )

met_cancer_freq = merge(cancer_freq, met_cancer_mut)
ex14_cancer_freq = merge(cancer_freq, ex14_cancer_mut)

plot_met_cancer_freq <- ggplot(met_cancer_freq, aes(x=IL3_withdrawal_score, y=cBioPortal_frq, label=hgvs))+
  geom_point()+
  xlab("Activity Score")+
  ylab("cBioPortal Mutation Frequency")+
  ggtitle("METdEx14")+
  xlim(-8,2)+
  geom_vline(xintercept=met_syn_mean, linetype="dotted",color="red")+
  theme_bw()

plot_ex14_cancer_freq <- ggplot(ex14_cancer_freq, aes(x=IL3_withdrawal_score, y=cBioPortal_frq, label=hgvs))+
  geom_point()+
  xlab("Activity Score")+
  ylab("cBioPortal Mutation Frequency")+
  ggtitle("MET+Ex14")+
  xlim(-8,2)+
  geom_vline(xintercept=ex14_syn_mean, linetype="dotted",color="red")+
  theme_bw()
plot(plot_met_cancer_freq)

plot_grid(plot_ex14_cancer_freq,plot_met_cancer_freq, ncol=1, nrow=2)


#####################################################################################################
###---------------- plot MET resistance mutations as distributions ----------------------------------
#####################################################################################################

### specific mutations 

met_resis_mut = met_scores_  %>% filter((hgvs %in% c("p.(D170N)","p.(D170H)","p.(Y172H)",
                                                    "p.(D170A)","p.(D170S)","p.(L137F)","p.(F142I)",
                                                    "p.(H36Y)","p.(L137V)","p.(Y172S)","p.(Y172C)",
                                                    "p.(V34L)","p.(G105R)","p.(F142L)")))

ex14_resis_mut = ex14_scores_  %>% filter((hgvs %in% c("p.(D170N)","p.(D170H)","p.(Y172H)",
                                                    "p.(D170A)","p.(D170S)","p.(L137F)","p.(F142I)",
                                                    "p.(H36Y)","p.(L137V)","p.(Y172S)","p.(Y172C)",
                                                    "p.(V34L)","p.(G105R)","p.(F142L)")))

met_resis_mut_df <-data.frame(hgvs=met_resis_mut$hgvs,
                              met_score=met_resis_mut$IL3_withdrawal_score,
                              mutation_type=met_resis_mut$mutation_type,
                              pos=met_resis_mut$pos)

ex14_resis_mut_df <-data.frame(hgvs=ex14_resis_mut$hgvs,
                               ex14_score=ex14_resis_mut$IL3_withdrawal_score,       
                               mutation_type=ex14_resis_mut$mutation_type, 
                               pos=ex14_resis_mut$pos)

merged_resis_mut_datasets <-merge(met_resis_mut_df,ex14_resis_mut_df)

plot <- ggplot(merged_resis_mut_datasets) +
  geom_density(aes(x = met_score)) +
  xlab("Enrich2 score, MET KD TKI reistance positions") + 
  ylab("Density") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_color_manual(values=cbp1)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")

plot(plot)

plot <- ggplot(merged_resis_mut_datasets) +
  geom_density(aes(x = ex14_score)) +
  xlab("Enrich2 score, Ex14 KD TKI reistance positions") + 
  ylab("Density") + 
  theme_classic() + 
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_color_manual(values=cbp1)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")

plot(plot)

#####################################################################################################
###---------------- plot MET resistance Positions as distributions ----------------------------------
#####################################################################################################
                                       
met_resis_pos = met_scores_  %>% filter((pos %in% c(1230,1200,1228,1092,1195,1163,1094,1164,1163)))

ex14_resis_pos = ex14_scores_ %>% filter((pos %in% c(1230,1200,1228,1092,1195,1163,1094,1164,1163)))

met_resis_pos_df<-data.frame(hgvs=met_resis_pos$hgvs,
                              met_score=met_resis_pos$IL3_withdrawal_score,
                              mutation_type=met_resis_pos$mutation_type,
                              pos=met_resis_pos$pos)

ex14_resis_pos_df<-data.frame(hgvs=ex14_resis_pos$hgvs,
                               ex14_score=ex14_resis_pos$IL3_withdrawal_score,       
                               mutation_type=ex14_resis_pos$mutation_type, 
                               pos=ex14_resis_pos$pos)

merged_resis_pos_datasets<-merge(met_resis_pos_df,ex14_resis_pos_df)

plot_met_resistance_pos <- ggplot(met_resis_pos) +
  geom_histogram(aes(x = IL3_withdrawal_score),color="black",fill="darkseagreen3", alpha=0.5) +
  xlab("Acivity Score") + 
  ylab("Count") + 
  ylim(0,30)+
  xlim(-10,3)+
  geom_vline(xintercept=met_syn_mean, linetype="dotted",color="red")+
  geom_text(aes(x=met_syn_mean-4, label="mean = -0.43", y=125), colour="red", angle=0)+
  ggtitle("METdEx14 resistance positions")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_met_resistance_pos)

plot_ex14_resistance_pos <- ggplot(ex14_resis_pos) +
  geom_histogram(aes(x = IL3_withdrawal_score),color="black",fill="darkseagreen3", alpha=0.5) +
  xlab("Activity score") + 
  ylim(0,30)+
  xlim(-10,3)+
  geom_vline(xintercept=ex14_syn_mean, linetype="dotted",color="red")+
  geom_text(aes(x=ex14_syn_mean-4, label="mean = -0.14", y=125), colour="red", angle=0)+
  ggtitle("MET+Ex14 resistance positions")+
  ylab("Count") +
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_ex14_resistance_pos)


### this code block filters out mutations at KNOWN RESISTANCE positions that are 2 SD above WT mean and have low error 
met_resistance_2SD <- subset(met_resis_pos %>% 
                           filter(
                             (abs(IL3_withdrawal_score - mean(met_syn$IL3_withdrawal_score)) >= 2*sd(met_syn$IL3_withdrawal_score)) 
                             & (IL3_withdrawal_SE <= met_wt_SE_mean))
                           )

plot_met_resis_scatter <- ggplot(met_resistance_2SD %>% filter(IL3_withdrawal_score>0),aes(x=pos,y=IL3_withdrawal_score, label=hgvs)) +
  geom_point()+
  ggtitle("METdEx14 GOF resistance pos")+
  geom_text(hjust=0, vjust=0)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())


ex14_resistance_2SD <- subset(ex14_resis_pos %>% 
                            filter(
                              (abs(IL3_withdrawal_score - mean(ex14_syn$IL3_withdrawal_score)) >= 2*sd(ex14_syn$IL3_withdrawal_score)) 
                              & (IL3_withdrawal_SE <= ex14_wt_SE_mean))
                          )


plot_ex14_resis_scatter <- ggplot(ex14_resistance_2SD %>% filter(IL3_withdrawal_score>0),aes(x=pos,y=IL3_withdrawal_score, label=hgvs)) +
  geom_point()+
  ggtitle("MET+Ex14 GOF resisatnce pos")+
  geom_text(hjust=0, vjust=0)+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

plot_grid(plot_ex14_resistance_pos,plot_ex14_resis_scatter,plot_met_resistance_pos,plot_met_resis_scatter,ncol=2,nrow=2)

plot_grid(plot_ex14_resistance_pos,plot_ecan, plot_ecan_scatter,
          plot_met_resistance_pos, plot_mcan, plot_mcan_scatter,
          ncol=3, nrow=2)


##### this code block plots the KNOWN resistance mutation distributions compared to total distributions 

cbp3 <- c("#999999","#00A693")
data5 = met_resis_mut$IL3_withdrawal_score
data6 = met_missense$IL3_withdrawal_score
dd = data.frame(x = c(data5,data6), 
               type=rep(c("Known resistance mutations", "All Mutations"), 
                        c(length(data5), length(data6))))

plot_met_resis_1 <- ggplot(dd, aes(x, color=type, fill=type)) +
  scale_colour_manual(values=cbp3)+
  scale_fill_manual(values=cbp3)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.6,
                 mapping = aes(y = stat(ncount)))+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 resistance mutations")+
  xlab("Activity Score")+
  ylab("Normalizaed Count")+
  xlim(-10,3)+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")


data7 = ex14_resis_mut$IL3_withdrawal_score
data8 = ex14_missense$IL3_withdrawal_score
cc = data.frame(x = c(data7,data8), 
               type=rep(c("Known resistance mutations", "All mutations"), 
                        c(length(data7), length(data8))))

plot_ex14_resis_1 <- ggplot(cc, aes(x, color=type, fill=type)) +
  scale_colour_manual(values=cbp3)+
  scale_fill_manual(values=cbp3)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.6,
                 mapping = aes(y = stat(ncount)))+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 resistance mutations")+
  xlab("Activity Score")+
  ylab("Normalizaed Count")+
  xlim(-10,3)+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")

# plots reistance mutaiton distibutions 
plot_grid(plot_ex14_resis_1,plot_met_resis_1, ncol=2,nrow=1)


#####################################################################################################
###--------------------------- Map Resistance Mutations 3D Structure  -------------------------------
#####################################################################################################

met_resis_pos_df_2 <- met_resis_pos_df %>% filter(met_score>0 & met_score > met_syn_mean)
met_Resistance <- data.frame(pos = met_resis_pos_df_2$pos,
                             score=met_resis_pos_df_2$met_score)
mean_met_Resistance <- met_Resistance %>% group_by(pos) %>% summarise(resis_avg = mean(score))
met_avg_Resistance_3R7O <- read.pdb("3R7O")x
x = map_scores_pdb(met_avg_Resistance_3R7O ,mean_met_Resistance, "resis_avg")
write.pdb(x, file="met_avg_Resistance_3R7O.pdb")

ex14_resis_pos_df_2 <- ex14_resis_pos_df %>% filter(ex14_score>0 & ex14_score > ex14_syn_mean)
ex14_Resistance <- data.frame(pos = ex14_resis_pos_df_2$pos,
                             score=ex14_resis_pos_df_2$ex14_score)
mean_ex14_Resistance <- ex14_Resistance %>% group_by(pos) %>% summarise(resis_avg = mean(score))
ex14_avg_Resistance_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_avg_Resistance_3R7O ,mean_ex14_Resistance, "resis_avg")
write.pdb(x, file="ex14_avg_Resistance_3R7O.pdb")

#####################################################################################################
###--------------------------- Resistance spatial differences  --------------------------------------
#####################################################################################################

# one of our key expectations is that resistance mutations that appear GOF in our screen are
# in relation to the catalytic site...not necessarily allosteric 
# to this I'm going to map the scores of the resitance mutations on structures assocaited w/ resistance 
# I'm going to do this for Crizotinib since there is a strcutrue...there are not strucures for all 

# crizotinib is 2WGJ 

ex14_crizotinib = ex14_scores_  %>% filter((hgvs %in% c("p.(D170N)",
                                                       "p.(D170H)",
                                                       "p.(D170S)",
                                                       "p.(D170A)",
                                                       "p.(Y172H)",
                                                       "p.(Y172S)",
                                                       "p.(Y172C)",
                                                       "p.(G105R)",
                                                       "p.(F142I)"
                                                       )))

ex14_crizotinib_avg_scores <- ex14_crizotinib %>% group_by(pos) %>% summarise(ex14_crizotinib_avg = mean(IL3_withdrawal_score, na.rm=TRUE))

pdb_ex14_crizotinib <- read.pdb("2WGJ") # open PDB
x = map_scores_pdb(pdb_ex14_crizotinib, ex14_crizotinib_avg_scores  , "ex14_crizotinib_avg")
write.pdb(x, file="ex14_crizotinib_avg_scores_2WGJ.pdb")
# spectrum b, firebrick_white_, minimum =-2, maximum=0


# merestinib 4EEV
ex14_merestinib = ex14_scores_  %>% filter((hgvs %in% c("p.(L137F)",
                                                        "p.(L137V)",
                                                        "p.(F142I)"
                                                        )))

ex14_merestinib_avg_scores <- ex14_merestinib %>% group_by(pos) %>% summarise(ex14_merestinib_avg = mean(IL3_withdrawal_score, na.rm=TRUE))

pdb_ex14_merestinib <- read.pdb("4EEV") # open PDB
x = map_scores_pdb(pdb_ex14_merestinib, ex14_merestinib_avg_scores  , "ex14_merestinib_avg")
write.pdb(x, file="ex14_merestinib_avg_scores_4EEV.pdb")
# spectrum b, firebrick_white_, minimum =-2, maximum=0



#####################################################################################################
###-------------- Hamming distance calculations for cancer and resistance mutations  ----------------
#####################################################################################################

# the goal of this code is to show that GOF mutations that exhibity higher activity scores than 
# assocaited, known cancer or resistance mutations are due to the experimental access to the genetic code 
# in other ways, mutations that are MORE GOF than known cancer/ resistance mutations are >2 bp away in the condon 

met_cancer_GOF = met_cancer_2SD %>% filter(IL3_withdrawal_score>0)

#met_resistance_2SD
#ex14_resistance_2SD

high_usage_codon_table <- data.frame (aa  = c('F','L','Y','STOP','H','Q','I','M','N','K','V','D','E','S','C','W','P','R','T','A','G'),
                  codon= c('TTC','CTG','TAC','TGA','CAC','CAG','ATC','ATG','AAC','AAG','GTG','GAC','GAG','AGC','TGC','TGG','CCC','CGG','ACC','GCC','GGC')
                  )

# creates base table of mutations that score better than known cancer mutations 
met_screen_cancer_GOF <- data.frame( met_screen_GOF_hgvs = met_cancer_GOF$hgvs,
                                     met_screen_GOF_score = met_cancer_GOF$IL3_withdrawal_score,
                                     met_screen_GOF_pos = met_cancer_GOF$pos,
                                     met_screen_GOF_variant_aa = met_cancer_GOF$variants
                                     )

met_screen_cancer_GOF <- met_screen_cancer_GOF  %>% filter((!met_screen_GOF_hgvs %in% c("p.(A193T)","p.(C33Y)","p.(C250S)","p.(D41G)","p.(D122Y)","p.(D228G)",
                                                                                                      "p.(E256G)","p.(E267K)","p.(F66L)","p.(H116N)","p.(H180L)",
                                                                                                      "p.(K141T)","p.(L82F)","p.(L128F)","p.(L128I)","p.(L128R)","p.(L147V)",
                                                                                                      "p.(L254V)","p.(M102I)","p.(M134I)","p.(N42D)","p.(N80K)","p.(P2S)",
                                                                                                      "p.(P15S)","p.(P95Q)","p.(P253L)","p.(Q65R)","p.(Q200K)","p.(R90Q)",
                                                                                                      "p.(R221I)","p.(R221K)","p.(R269C)","p.(S178R)","p.(T38S)","p.(V5A)",
                                                                                                      "p.(V12E)","p.(V12G)","p.(V12M)","p.(V63A)","p.(V229I)","p.(Y177H)","p.(M192T)", 
                                                                                                      "p.(M192I)",  "p.(Y172S)", "p.(Y172C)","p.(Y172N)", "p.(Y172H)",
                                                                                                      "p.(F142I)", "p.(H36Y)","p.(D170N)", "p.(D170A)","p.(D170S)","p.(D170H)",
                                                                                                      "p.(V34I)", "p.(G105R)","p.(L137V)","p.(L137F)")))


# creates base table of known cancer mutations 
met_known_cancer_codons <- data.frame( met_known_cancer_hgvs = met_cancer_mut$hgvs,
                                     met_known_cancer_score = met_cancer_mut$IL3_withdrawal_score,
                                     met_known_cancer_pos = met_cancer_mut$pos,
                                     met_known_cancer_variant_aa = met_cancer_mut$variants
)

#### creates a table with the mutations that are MORE GOF than the detected cancer 
met_screen_cancer_GOF['wt_aa'] <- NA
for (i in met_screen_cancer_GOF){
  x = substr(met_screen_cancer_GOF$met_screen_GOF_hgvs, 4,4) # this is the WT codon 
  met_screen_cancer_GOF['wt_aa'] <- x
}

met_screen_cancer_GOF = met_screen_cancer_GOF %>%
  left_join(high_usage_codon_table, by = c(wt_aa= "aa")) %>%
  left_join(high_usage_codon_table, by = c(met_screen_GOF_variant_aa= "aa")) %>%
  left_join(met_known_cancer_codons, by = c(met_screen_GOF_pos= "met_known_cancer_pos"))  

met_screen_cancer_GOF = met_screen_cancer_GOF %>%
  left_join(high_usage_codon_table, by = c(met_known_cancer_variant_aa= "aa"))


colnames(met_screen_cancer_GOF)[6] ="wt_codon"
colnames(met_screen_cancer_GOF)[7] ="GOF_variant_codon"
colnames(met_screen_cancer_GOF)[11] ="cancer_variant_codon"


### function to compute hamming distance 
# distance between cancer_variant_codon and GOF_variant_codon
# vs distance between cancer_variant_codon and wt_codon 

x <- met_screen_cancer_GOF$GOF_variant_codon
y <- met_screen_cancer_GOF$cancer_variant_codon
z <- met_screen_cancer_GOF$wt_codon

## hamming distance
p <- stringdist(z,x,method="hamming")  # wt to GOF
q <- stringdist(z,y,method="hamming")  # wt to cancer

met_screen_cancer_GOF['wt_vs_GOF_dist'] <- p
met_screen_cancer_GOF['cancer_vs_wt_dist'] <- q


met_screen_cancer_GOF = met_screen_cancer_GOF %>% filter(wt_vs_GOF_dist != 0)
#### split the above data frame to eliminate double counting 
met_known_cancer = met_screen_cancer_GOF[!duplicated(met_screen_cancer_GOF$met_known_cancer_hgvs),]

### hamming distribution plot 
plot_met_hamming_1 <- ggplot()+ 
  geom_histogram(aes(x=met_known_cancer$cancer_vs_wt_dist),alpha=0.7,binwidth=0.5,fill="#00B2CA",color="#00B2CA")+
  theme_classic()+
  xlab("Hamming distance")+
  ggtitle("METdEx14")+
  ylim(0,5)+
  xlim(0,3.5)+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot_met_hamming_2 <- ggplot()+ 
  geom_histogram(aes(x=met_known_cancer$wt_vs_GOF_dist),alpha=0.7,binwidth=0.5,fill="#999999",color="#999999")+
  theme_classic()+
  xlab("Hamming distance")+
  ggtitle("METdEx14")+
  ylim(0,5)+
  xlim(0,3.5)+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot_grid(plot_met_hamming_1,plot_met_hamming_2,ncol=2,nrow=1)

  

#####################################################################
# same but now for Ex14

ex14_cancer_GOF = ex14_cancer_2SD %>% filter(IL3_withdrawal_score>0)

#ex14_resistance_2SD
#ex14_resistance_2SD

high_usage_codon_table <- data.frame (aa  = c('F','L','Y','STOP','H','Q','I','M','N','K','V','D','E','S','C','W','P','R','T','A','G'),
                                      codon= c('TTC','CTG','TAC','TGA','CAC','CAG','ATC','ATG','AAC','AAG','GTG','GAC','GAG','AGC','TGC','TGG','CCC','CGG','ACC','GCC','GGC')
)

# creates base table of mutations that score better than known cancer mutations 
ex14_screen_cancer_GOF <- data.frame( ex14_screen_GOF_hgvs = ex14_cancer_GOF$hgvs,
                                     ex14_screen_GOF_score = ex14_cancer_GOF$IL3_withdrawal_score,
                                     ex14_screen_GOF_pos = ex14_cancer_GOF$pos,
                                     ex14_screen_GOF_variant_aa = ex14_cancer_GOF$variants
)

ex14_screen_cancer_GOF <- ex14_screen_cancer_GOF  %>% filter((!ex14_screen_GOF_hgvs %in% c("p.(A193T)","p.(C33Y)","p.(C250S)","p.(D41G)","p.(D122Y)","p.(D228G)",
                                                                                        "p.(E256G)","p.(E267K)","p.(F66L)","p.(H116N)","p.(H180L)",
                                                                                        "p.(K141T)","p.(L82F)","p.(L128F)","p.(L128I)","p.(L128R)","p.(L147V)",
                                                                                        "p.(L254V)","p.(M102I)","p.(M134I)","p.(N42D)","p.(N80K)","p.(P2S)",
                                                                                        "p.(P15S)","p.(P95Q)","p.(P253L)","p.(Q65R)","p.(Q200K)","p.(R90Q)",
                                                                                        "p.(R221I)","p.(R221K)","p.(R269C)","p.(S178R)","p.(T38S)","p.(V5A)",
                                                                                        "p.(V12E)","p.(V12G)","p.(V12M)","p.(V63A)","p.(V229I)","p.(Y177H)","p.(M192T)", 
                                                                                        "p.(M192I)",  "p.(Y172S)", "p.(Y172C)","p.(Y172N)", "p.(Y172H)",
                                                                                        "p.(F142I)", "p.(H36Y)","p.(D170N)", "p.(D170A)","p.(D170S)","p.(D170H)",
                                                                                        "p.(V34I)", "p.(G105R)","p.(L137V)","p.(L137F)")))

# creates base table of known cancer mutations 
ex14_known_cancer_codons <- data.frame( ex14_known_cancer_hgvs = ex14_cancer_mut$hgvs,
                                       ex14_known_cancer_score = ex14_cancer_mut$IL3_withdrawal_score,
                                       ex14_known_cancer_pos = ex14_cancer_mut$pos,
                                       ex14_known_cancer_variant_aa = ex14_cancer_mut$variants
)

#### creates a table with the mutations that are MORE GOF than the known cancer mutation 
ex14_screen_cancer_GOF['wt_aa'] <- NA
for (i in ex14_screen_cancer_GOF){
  x = substr(ex14_screen_cancer_GOF$ex14_screen_GOF_hgvs, 4,4) # this is the WT codon 
  ex14_screen_cancer_GOF['wt_aa'] <- x
}

ex14_screen_cancer_GOF = ex14_screen_cancer_GOF %>%
  left_join(high_usage_codon_table, by = c(wt_aa= "aa")) %>%
  left_join(high_usage_codon_table, by = c(ex14_screen_GOF_variant_aa= "aa")) %>%
  left_join(ex14_known_cancer_codons, by = c(ex14_screen_GOF_pos= "ex14_known_cancer_pos"))  

ex14_screen_cancer_GOF = ex14_screen_cancer_GOF %>%
  left_join(high_usage_codon_table, by = c(ex14_known_cancer_variant_aa= "aa"))


colnames(ex14_screen_cancer_GOF)[6] ="wt_codon"
colnames(ex14_screen_cancer_GOF)[7] ="GOF_variant_codon"
colnames(ex14_screen_cancer_GOF)[11] ="cancer_variant_codon"


### function to compute hamming distance 
# distance between cancer_variant_codon and GOF_variant_codon
# vs distance between cancer_variant_codon and wt_codon 

x <- ex14_screen_cancer_GOF$GOF_variant_codon
y <- ex14_screen_cancer_GOF$cancer_variant_codon
z <- ex14_screen_cancer_GOF$wt_codon

## hamming distance
p <- stringdist(z,x,method="hamming")  # wt to GOF
q <- stringdist(z,y,method="hamming")  # wt to cancer

ex14_screen_cancer_GOF['wt_vs_GOF_dist'] <- p
ex14_screen_cancer_GOF['cancer_vs_wt_dist'] <- q


ex14_screen_cancer_GOF = ex14_screen_cancer_GOF %>% filter(wt_vs_GOF_dist != 0)
#### split the above data frame to eliminate double counting 
ex14_known_cancer = ex14_screen_cancer_GOF[!duplicated(ex14_screen_cancer_GOF$ex14_known_cancer_hgvs),]

### hamming distribution plot 
plot_ex14_hamming_1 <- ggplot()+ 
  geom_histogram(aes(x=ex14_known_cancer$cancer_vs_wt_dist),alpha=0.7,binwidth=0.5,fill="#00B2CA",color="#00B2CA")+
  theme_classic()+
  xlab("Hamming distance")+
  ggtitle("MET+Ex14")+
  ylim(0,13)+
  xlim(0,3.5)+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot_ex14_hamming_2 <- ggplot()+ 
  geom_histogram(aes(x=ex14_known_cancer$wt_vs_GOF_dist),alpha=0.7,binwidth=0.5,fill="#999999",color="#999999")+
  theme_classic()+
  xlab("Hamming distance")+
  ggtitle("MET+Ex14")+
  ylim(0,13)+
  xlim(0,3.5)+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot_grid(plot_ex14_hamming_1 ,plot_ex14_hamming_2, ncol=2,nrow=1)

#####################################################################################################
###---------------- Map Hamming Distance Residues on the structures  -----------------------------
#####################################################################################################

### real data 
ex14_avg_GOF_hamming <- data.frame(pos = ex14_known_cancer$ex14_screen_GOF_pos,HD = ex14_known_cancer$wt_vs_GOF_dist)
ex14_avg_GOF_hamming_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_avg_GOF_hamming_3R7O, ex14_avg_GOF_hamming, "HD")
write.pdb(x, file="ex14_avg_GOF_hamming_3R7O.pdb")

met_avg_GOF_hamming <- data.frame(pos = met_known_cancer$met_screen_GOF_pos,HD = met_known_cancer$wt_vs_GOF_dist)
met_avg_GOF_hamming_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_avg_GOF_hamming_3R7O, met_avg_GOF_hamming, "HD")
write.pdb(x, file="met_avg_GOF_hamming_3R7O.pdb")

#### Creates Hamming distance legend, where ribbon thickness = HD 

HD_legend0 <- data.frame(pos = c(1200,1293), HD = c(3,0))
HD_legend1 <- data.frame(pos = c(1200,1293), HD = c(3,1))
HD_legend2 <- data.frame(pos = c(1200,1293), HD = c(3,2))
HD_legend3 <- data.frame(pos = c(1200,1293), HD = c(1,3))
HD_legend <- read.pdb("3R7O")
x0 = map_scores_pdb(HD_legend, HD_legend0 , "HD")
x1 = map_scores_pdb(HD_legend, HD_legend1 , "HD")
x2 = map_scores_pdb(HD_legend, HD_legend2 , "HD")
x3 = map_scores_pdb(HD_legend, HD_legend3 , "HD")
write.pdb(x0, file="HD_legend0_3R7O.pdb")
write.pdb(x1, file="HD_legend1_3R7O.pdb")
write.pdb(x2, file="HD_legend2_3R7O.pdb")
write.pdb(x3, file="HD_legend3_3R7O.pdb")


#####################################################################################################
###---------------- Mutations Detected, Validated, and GOF  -----------------------------------------
#####################################################################################################

### double check positions
met_detected <- met_scores_  %>% filter((hgvs %in% c("p.(A193T)","p.(C33Y)","p.(C250S)","p.(D41G)","p.(D122Y)","p.(D228G)",
                                                     "p.(E256G)","p.(E267K)","p.(F66L)","p.(H116N)","p.(H180L)",
                                                     "p.(K141T)","p.(L82F)","p.(L128F)","p.(L128I)","p.(L128R)","p.(L147V)",
                                                     "p.(L254V)","p.(M102I)","p.(M134I)","p.(N42D)","p.(N80K)","p.(P2S)",
                                                     "p.(P15S)","p.(P95Q)","p.(P253L)","p.(Q65R)","p.(Q200K)","p.(R90Q)",
                                                     "p.(R221I)","p.(R221K)","p.(R269C)","p.(S178R)","p.(T38S)","p.(V5A)",
                                                     "p.(V12E)","p.(V12G)","p.(V12M)","p.(V63A)","p.(V229I)","p.(Y177H)")))

met_validated <- met_scores_  %>% filter((hgvs %in% c("p.(M192T)", "p.(M192I)",  "p.(Y172S)", "p.(Y172C)","p.(Y172N)", "p.(Y172H)",
                                                        "p.(F142I)", "p.(H36Y)","p.(D170N)", "p.(D170A)","p.(D170S)","p.(D170H)",
                                                        "p.(V34I)", "p.(G105R)","p.(L137V)","p.(L137F)")))


met_GOF_better_than_detected <- data.frame( met_screen_GOF_hgvs = met_cancer_GOF$hgvs,
                                             met_screen_GOF_score = met_cancer_GOF$IL3_withdrawal_score,
                                             met_screen_GOF_pos = met_cancer_GOF$pos,
                                             met_screen_GOF_variant_aa = met_cancer_GOF$variants)

met_GOF_better_than_detected <- met_GOF_better_than_detected  %>% filter((!met_screen_GOF_hgvs %in% c("p.(A193T)","p.(C33Y)","p.(C250S)","p.(D41G)","p.(D122Y)","p.(D228G)",
                                                                                      "p.(E256G)","p.(E267K)","p.(F66L)","p.(H116N)","p.(H180L)",
                                                                                      "p.(K141T)","p.(L82F)","p.(L128F)","p.(L128I)","p.(L128R)","p.(L147V)",
                                                                                      "p.(L254V)","p.(M102I)","p.(M134I)","p.(N42D)","p.(N80K)","p.(P2S)",
                                                                                      "p.(P15S)","p.(P95Q)","p.(P253L)","p.(Q65R)","p.(Q200K)","p.(R90Q)",
                                                                                      "p.(R221I)","p.(R221K)","p.(R269C)","p.(S178R)","p.(T38S)","p.(V5A)",
                                                                                      "p.(V12E)","p.(V12G)","p.(V12M)","p.(V63A)","p.(V229I)","p.(Y177H)","p.(M192T)", 
                                                                                      "p.(M192I)",  "p.(Y172S)", "p.(Y172C)","p.(Y172N)", "p.(Y172H)",
                                                                                      "p.(F142I)", "p.(H36Y)","p.(D170N)", "p.(D170A)","p.(D170S)","p.(D170H)",
                                                                                      "p.(V34I)", "p.(G105R)","p.(L137V)","p.(L137F)")))



######### plots GOF mutations that are stronger than validated cancer mutations 

met_detected_better_than_validated <- met_detected %>% filter(IL3_withdrawal_score>0.5)

plot_met_detected_better_than_validated <-ggplot()+
  geom_point(aes(x=met_detected_better_than_validated$IL3_withdrawal_score,
                 y=met_detected_better_than_validated$pos))+
  geom_text(aes(x=met_detected_better_than_validated$IL3_withdrawal_score,
                y=met_detected_better_than_validated$pos,
                label=met_detected_better_than_validated$hgvs),nudge_x = 0.15, nudge_y = 0.25)+
  xlab("Activity Score")+
  ylab("Position")+
  ggtitle("METdEx14")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot(plot_met_detected_better_than_validated)


### detected cancer mutations 
met_ded = met_detected$IL3_withdrawal_score # detected distribution 
met_miss = met_missense$IL3_withdrawal_score #missense distribution 
met_ded_miss = data.frame(x = c(met_ded,met_miss), 
                         type=rep(c("met_detected", "All mutations"), 
                                  c(length(met_ded ), length(met_miss))))
plot_met_detected<- ggplot(met_ded_miss,aes(x, color=type, fill=type))+ 
  scale_colour_manual(values=cbp1)+
  scale_fill_manual(values=cbp1)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.5,
                 mapping = aes(y = stat(ncount)))+
  #geom_point(aes(y=met_detected$pos,x=met_detected$IL3_withdrawal_score))+
  geom_vline(xintercept=met_syn_mean, linetype="dashed")+
  ylab("Normalized Count")+
  xlab("Activity Score")+
  xlim(-10,5)+
  ggtitle("met all detected")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot(plot_met_detected)

### validated mutations 
met_val = met_validated$IL3_withdrawal_score # validated distribution 
met_val_miss = data.frame(x = c(met_val,met_miss), 
                         type=rep(c("met_val", "All mutations"), 
                                  c(length(met_val ), length(met_miss))))

plot_met_cancer_validated <- ggplot(met_val_miss,aes(x, color=type, fill=type))+ 
  scale_colour_manual(values=cbp1)+
  scale_fill_manual(values=cbp1)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.5,
                 mapping = aes(y = stat(ncount)))+
  #geom_point(aes(y=met_validated$pos,x=met_validated$IL3_withdrawal_score))+
  geom_vline(xintercept=met_syn_mean, linetype="dashed")+
  ylab("Normalized Count")+
  xlab("Activity Score")+
  xlim(-10,5)+
  ggtitle("met validated cancer")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot(plot_met_detected)
plot(plot_met_cancer_validated )

plot_met_GOF_better_than_detected <- ggplot()+ 
  geom_point(aes(y=met_GOF_better_than_detected$met_screen_GOF_pos,x=met_GOF_better_than_detected$met_screen_GOF_score))+
  geom_text(aes(y=met_GOF_better_than_detected$met_screen_GOF_pos,
                x=met_GOF_better_than_detected$met_screen_GOF_score,
                label=met_GOF_better_than_detected$met_screen_GOF_hgvs))+
  ylab("Position")+
  xlab("Activity Score")+
  ggtitle("MetdEx14 GOF> all cancer-assocaited")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot(plot_met_GOF_better_than_detected )


######### plots validated vs not validated cancer mutations 
cbp2 <- c("#F79256","#00B2CA","#999999")

met_GOF_better_cancer = met_GOF_better_than_detected$met_screen_GOF_score
met_ded_val_GOF = data.frame(x = c(met_ded, met_val, met_GOF_better_cancer), 
                          type=rep(c("cancer detected", "cancer validated","GOF"), 
                                   c(length(met_ded), length(met_val), length(met_GOF_better_cancer))))
plot_met_validated_vs_detected <- ggplot(met_ded_val_GOF,aes(x, color=type, fill=type))+ 
  scale_colour_manual(values=cbp2)+
  scale_fill_manual(values=cbp2)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.6,
                 mapping = aes(y = stat(ncount)))+
  #geom_point(aes(y=met_validated$pos,x=met_validated$IL3_withdrawal_score))+
  geom_vline(xintercept=met_syn_mean, linetype="dashed")+
  ylab("Normalized Count")+
  xlab("Activity Score")+
  xlim(-10,5)+
  ggtitle("METdEx14")+
  theme_classic()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot(plot_met_validated_vs_detected)



plot_met_det_val_dot <- ggplot()+
  geom_point((aes(x=met_detected$IL3_withdrawal_score,y=met_detected$pos, color="detected")))+
  geom_point((aes(x=met_validated$IL3_withdrawal_score,y=met_validated$pos, color="validated")))+
  geom_vline(xintercept=met_syn_mean,linetype="dashed")
plot(plot_met_det_val_dot)
  
plot_grid(plot_met_detected,plot_met_cancer_validated,plot_met_GOF_better_than_detected,ncol=3,nrow=1)




### double check positions
ex14_detected <- ex14_scores_  %>% filter((hgvs %in% c("p.(A193T)","p.(C33Y)","p.(C250S)","p.(D41G)","p.(D122Y)","p.(D228G)",
                                                       "p.(E256G)","p.(E267K)","p.(F66L)","p.(H116N)","p.(H180L)",
                                                       "p.(K141T)","p.(L82F)","p.(L128F)","p.(L128I)","p.(L128R)","p.(L147V)",
                                                       "p.(L254V)","p.(M102I)","p.(M134I)","p.(N42D)","p.(N80K)","p.(P2S)",
                                                       "p.(P15S)","p.(P95Q)","p.(P253L)","p.(Q65R)","p.(Q200K)","p.(R90Q)",
                                                       "p.(R221I)","p.(R221K)","p.(R269C)","p.(S178R)","p.(T38S)","p.(V5A)",
                                                       "p.(V12E)","p.(V12G)","p.(V12M)","p.(V63A)","p.(V229I)","p.(Y177H)")))

ex14_validated <- ex14_scores_  %>% filter((hgvs %in% c("p.(M192T)", "p.(M192I)",  "p.(Y172S)", "p.(Y172C)","p.(Y172N)", "p.(Y172H)",
                                                      "p.(F142I)", "p.(H36Y)","p.(D170N)", "p.(D170A)","p.(D170S)","p.(D170H)",
                                                      "p.(V34I)", "p.(G105R)","p.(L137V)","p.(L137F)")))

ex14_GOF_better_than_detected <- data.frame( ex14_screen_GOF_hgvs = ex14_cancer_GOF$hgvs,
                                            ex14_screen_GOF_score = ex14_cancer_GOF$IL3_withdrawal_score,
                                            ex14_screen_GOF_pos = ex14_cancer_GOF$pos,
                                            ex14_screen_GOF_variant_aa = ex14_cancer_GOF$variants)

ex14_GOF_better_than_detected <- ex14_GOF_better_than_detected  %>% filter((!ex14_screen_GOF_hgvs %in% c("p.(A193T)","p.(C33Y)","p.(C250S)","p.(D41G)","p.(D122Y)","p.(D228G)",
                                                                                                      "p.(E256G)","p.(E267K)","p.(F66L)","p.(H116N)","p.(H180L)",
                                                                                                      "p.(K141T)","p.(L82F)","p.(L128F)","p.(L128I)","p.(L128R)","p.(L147V)",
                                                                                                      "p.(L254V)","p.(M102I)","p.(M134I)","p.(N42D)","p.(N80K)","p.(P2S)",
                                                                                                      "p.(P15S)","p.(P95Q)","p.(P253L)","p.(Q65R)","p.(Q200K)","p.(R90Q)",
                                                                                                      "p.(R221I)","p.(R221K)","p.(R269C)","p.(S178R)","p.(T38S)","p.(V5A)",
                                                                                                      "p.(V12E)","p.(V12G)","p.(V12M)","p.(V63A)","p.(V229I)","p.(Y177H)","p.(M192T)", 
                                                                                                      "p.(M192I)",  "p.(Y172S)", "p.(Y172C)","p.(Y172N)", "p.(Y172H)",
                                                                                                      "p.(F142I)", "p.(H36Y)","p.(D170N)", "p.(D170A)","p.(D170S)","p.(D170H)",
                                                                                                      "p.(V34I)", "p.(G105R)","p.(L137V)","p.(L137F)")))



ex_ded = ex14_detected$IL3_withdrawal_score
ex_miss = ex14_missense$IL3_withdrawal_score
ex_ded_miss = data.frame(x = c(ex_ded,ex_miss), 
                              type=rep(c("ex14_detected", "All mutations"), 
                                       c(length(ex_ded ), length(ex_miss))))
plot_ex14_detected<- ggplot(ex_ded_miss,aes(x, color=type, fill=type))+ 
  scale_colour_manual(values=cbp1)+
  scale_fill_manual(values=cbp1)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.5,
                 mapping = aes(y = stat(ncount)))+
  #geom_point(aes(y=ex14_detected$pos,x=ex14_detected$IL3_withdrawal_score))+
  geom_vline(xintercept=ex14_syn_mean, linetype="dashed")+
  ylab("Normalized Count")+
  xlab("Activity Score")+
  xlim(-10,5)+
  ggtitle("ex14 all detected")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot(plot_ex14_detected)


ex_val = ex14_validated$IL3_withdrawal_score
ex_val_miss = data.frame(x = c(ex_val,ex_miss), 
                         type=rep(c("ex14_val", "All mutations"), 
                                  c(length(ex_val ), length(ex_miss))))

plot_ex14_cancer_validated <- ggplot(ex_val_miss,aes(x, color=type, fill=type))+ 
  scale_colour_manual(values=cbp1)+
  scale_fill_manual(values=cbp1)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.5,
                 mapping = aes(y = stat(ncount)))+
  #geom_point(aes(y=ex14_validated$pos,x=ex14_validated$IL3_withdrawal_score))+
  geom_vline(xintercept=ex14_syn_mean, linetype="dashed")+
  ylab("Normalized Count")+
  xlab("Activity Score")+
  xlim(-10,5)+
  ggtitle("ex14 validated cancer")+
  theme_bw()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot(plot_ex14_detected)
plot(plot_ex14_cancer_validated )

plot_ex14_GOF_better_than_detected <- ggplot()+ 
  geom_point(aes(y=ex14_GOF_better_than_detected$ex14_screen_GOF_pos,x=ex14_GOF_better_than_detected$ex14_screen_GOF_score))+
  geom_text(aes(y=ex14_GOF_better_than_detected$ex14_screen_GOF_pos,
                x=ex14_GOF_better_than_detected$ex14_screen_GOF_score,
                label=ex14_GOF_better_than_detected$ex14_screen_GOF_hgvs))+
  ylab("Position")+
  xlab("Activity Score")+
  ggtitle("MET+Ex14 GOF> all cancer-assocaited")+
  theme_bw()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")
plot(plot_ex14_GOF_better_than_detected )

plot_grid(plot_ex14_cancer_validated,plot_ex14_detected,plot_ex14_GOF_better_than_detected,
          plot_met_cancer_validated,plot_met_detected,plot_met_GOF_better_than_detected,
          ncol=3,nrow=2)


######### plots validated vs not validated cancer mutations 

ex14_GOF_better_cancer = ex14_GOF_better_than_detected$ex14_screen_GOF_score
ex_ded_val_GOF = data.frame(x = c(ex_ded,ex_val,ex14_GOF_better_cancer), 
                         type=rep(c("cancer detected", "cancer validated","GOF"), 
                                  c(length(ex_ded), length(ex_val),length(ex14_GOF_better_cancer))))
plot_ex14_validated_vs_detected <- ggplot(ex_ded_val_GOF,aes(x, color=type, fill=type))+ 
  scale_colour_manual(values=cbp2)+
  scale_fill_manual(values=cbp2)+
  geom_histogram(bins = 30, position = "identity", alpha = 0.6,
                 mapping = aes(y = stat(ncount)))+
  #geom_point(aes(y=met_validated$pos,x=met_validated$IL3_withdrawal_score))+
  geom_vline(xintercept=ex14_syn_mean, linetype="dashed")+
  ylab("Normalized Count")+
  xlab("Activity Score")+
  xlim(-10,5)+
  ggtitle("MET")+  
  theme_classic2()+
  theme(text = element_text(size = 15),
        #panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="none")

plot_grid(plot_ex14_validated_vs_detected ,plot_met_validated_vs_detected,ncol=1,nrow=2)

plot_grid(plot_ex14_validated_vs_detected ,plot_ex14_hamming_1, plot_ex14_hamming_2,
          plot_met_validated_vs_detected, plot_met_hamming_1, plot_met_hamming_2,
          ncol=3,nrow=2)


#####################################################################################################
###------------------ heatmap of GAIN of Function MetdEx14 mutations -------------------------------
#####################################################################################################

order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')

met_wt_1 = str_split(substr(met_wt_sequence, 1, 100), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 101, 200), '')[[1]]
met_wt_3 = str_split(substr(met_wt_sequence, 201, 287), '')[[1]]

met_row1 = ggplot(data = met_scores_ %>% filter(pos %in% c(1059:1158)),
                  aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(0,5)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1158, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1158),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

met_row2 = ggplot(data = met_scores_ %>% filter(pos %in% c(1159:1258)), 
                  aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(0,5)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1159,1258, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1159,1258),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

met_row3 = ggplot(data = met_scores_ %>% filter(pos %in% c(1259:1345)), 
                  aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(0,5)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1259,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1259,1345),
                       labels = met_wt_3,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right") +
  labs(y = "Mutation", x = "Position")

met_DMS = ggarrange(met_row1, met_row2, met_row3,
                    nrow = 3, ncol = 1)

ggsave("met_GOF_heatmap.pdf", height = 7, width = 8.5, met_DMS)
ggsave("met_GOF_heatmap.png", height = 7, width = 8.5, met_DMS)


#####################################################################################################
###----------- heatmap of LOSS of FUNCTION MetdEx14 mutations -------  ----------------------
#####################################################################################################

order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')

met_wt_1 = str_split(substr(met_wt_sequence, 1, 100), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 101, 200), '')[[1]]
met_wt_3 = str_split(substr(met_wt_sequence, 201, 287), '')[[1]]

met_row1 = ggplot(data = met_scores_ %>% filter(pos %in% c(1059:1158)),
                  aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-15,0)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1158, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1158),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

met_row2 = ggplot(data = met_scores_ %>% filter(pos %in% c(1159:1258)), 
                  aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-15,0)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1159,1258, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1159,1258),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

met_row3 = ggplot(data = met_scores_ %>% filter(pos %in% c(1259:1345)), 
                  aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-15,0)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1259,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1259,1345),
                       labels = met_wt_3,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right") +
  labs(y = "Mutation", x = "Position")

met_DMS = ggarrange(met_row1, met_row2, met_row3,
                    nrow = 3, ncol = 1)

ggsave("met_LOF_heatmap.pdf", height = 7, width = 8.5, met_DMS)
ggsave("met_LOF_heatmap.png", height = 7, width = 8.5, met_DMS)

#####################################################################################################
###----------- heatmap of GAIN of Function Met+Ex14 mutations ------- ----------------------
#####################################################################################################

order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')

ex14_wt_1 = str_split(substr(ex14_wt_sequence, 1, 100), '')[[1]]
ex14_wt_2 = str_split(substr(ex14_wt_sequence, 101, 200), '')[[1]]
ex14_wt_3 = str_split(substr(ex14_wt_sequence, 201, 287), '')[[1]]

ex14_row1 = ggplot(data = ex14_scores_ %>% filter(pos %in% c(1059:1158)),
                   aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(1,4)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1158, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1158),
                       labels = ex14_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

ex14_row2 = ggplot(data = ex14_scores_ %>% filter(pos %in% c(1159:1258)), 
                   aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(1,4)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1159,1258, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1159,1258),
                       labels = ex14_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

ex14_row3 = ggplot(data = ex14_scores_ %>% filter(pos %in% c(1259:1345)), 
                   aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(1,4)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1259,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1259,1345),
                       labels = ex14_wt_3,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right"
  ) +
  labs(y = "Mutation", x = "Position")

ex14_DMS = ggarrange(ex14_row1, ex14_row2, ex14_row3,
                     nrow = 3, ncol = 1)

ggsave("Ex14_GOF_heatmap.pdf", height = 7, width = 8.5, met_DMS)
ggsave("Ex14_GOF_heatmap.png", height = 7, width = 8.5, met_DMS)


###----------- heatmap of LOSS of FUNCTION Met+Ex14 mutations -------  ----------------------

                                       

order <- c('X','H', 'K', 'R', 'D', 'E',
           'C','M','N','Q','S','T','A','I','L','V','W','F',
           'Y','G','P')

variant_names <- c('Stop','H', 'K', 'R', 'D', 'E',
                   'C','M','N','Q','S','T','A','I','L','V','W','F',
                   'Y','G','P')

ex14_wt_1 = str_split(substr(ex14_wt_sequence, 1, 100), '')[[1]]
ex14_wt_2 = str_split(substr(ex14_wt_sequence, 101, 200), '')[[1]]
ex14_wt_3 = str_split(substr(ex14_wt_sequence, 201, 287), '')[[1]]

ex14_row1 = ggplot(data = ex14_scores_ %>% filter(pos %in% c(1059:1158)),
                   aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,0)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1158, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1158),
                       labels = ex14_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

ex14_row2 = ggplot(data = ex14_scores_ %>% filter(pos %in% c(1159:1258)), 
                   aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,0)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1159,1258, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1159,1258),
                       labels = ex14_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

ex14_row3 = ggplot(data = ex14_scores_ %>% filter(pos %in% c(1259:1345)), 
                   aes(x = pos, y = factor(variants, level = order, labels = variant_names), fill = IL3_withdrawal_score )) +
  geom_tile(aes(color = as.factor(is.wt)), size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,0)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1259,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1259,1345),
                       labels = ex14_wt_3,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="right"
  ) +
  labs(y = "Mutation", x = "Position")

ex14_DMS = ggarrange(ex14_row1, ex14_row2, ex14_row3,
                     nrow = 3, ncol = 1)

ggsave("Ex14_LOF_heatmap.pdf", height = 7, width = 8.5, met_DMS)
ggsave("Ex14_LOF_heatmap.png", height = 7, width = 8.5, met_DMS)



###----------- bin pos and average positional fitness scores for GOF and LOF  ----------------------

# IL3 withdrawal
GOF_df<-subset(merged_datasets, met_score > 0 & ex14_score > 0) # gain of function mutations 
LOF_df<-subset(merged_datasets, met_score < 0 & ex14_score < 0) # loss of function mutations

#### averages for GOF and LOF 
GOF_avg_df <- GOF_df %>% group_by(pos) %>% summarise(GOF_met = mean(met_score, na.rm=TRUE),
                                                     GOF_ex14 = mean(ex14_score, na.rm=TRUE))
LOF_avg_df <- LOF_df %>% group_by(pos) %>% summarise(LOF_met = mean(met_score, na.rm=TRUE),
                                                     LOF_ex14 = mean(ex14_score, na.rm=TRUE))

### max scores for GOF and LOF 
GOF_max_df <- GOF_df %>% group_by(pos) %>% summarise(max_GOF_met = max(met_score, na.rm=TRUE),
                                                     max_GOF_ex14 = max(ex14_score, na.rm=TRUE))
LOF_min_df <- LOF_df %>% group_by(pos) %>% summarise(min_LOF_met = min(met_score, na.rm=TRUE),
                                                     min_LOF_ex14 = min(ex14_score, na.rm=TRUE))

#####################################################################################################
##------------------------------------ Physiochemical mapping ---------------------------------------
#####################################################################################################


# hydrophobic: V,L,I,M,P,G,A,W,F
# aromatic : F,Y,W 
# acidic : D,E 
# basic : H,K,R
# amide : N,Q 
# nucleophilic: S,T,C
# polar, uncharged: S,T,Y,N,Q, C
# small: G, A
# structure breaking: P,G
# aliphatic: I,L,V,A,G
# sulfer containing: M, C


###------ global physiochem mapping ---------

####################### Hydrophobic 
hydrophobic = merged_datasets %>% filter(variants %in% c("V","L","I","M","P","G","A","W","F"))

plot3 <- ggplot(data=hydrophobic, aes(x=met_score, y=ex14_score, color = mutation_type)) + 
  geom_point(shape=1,aes(color = mutation_type)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("Hydrophobic")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 <- ggplot(hydrophobic, aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens4 <- ggplot(hydrophobic, aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))


############# aromatic
aromatic = merged_datasets %>% filter(variants %in% c("Y","W","F"))

plot3 <- ggplot(data=aromatic, aes(x=met_score, y=ex14_score, color = mutation_type)) + 
  geom_point(shape=1,aes(color = mutation_type)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("Aromatic")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 <- ggplot(aromatic, aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens4 <- ggplot(aromatic, aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

############# basic
basic = merged_datasets %>% filter(variants %in% c( "H","K","R"))

plot3 <- ggplot(data=basic, aes(x=met_score, y=ex14_score, color = mutation_type)) + 
  geom_point(shape=1,aes(color = mutation_type)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("Basic")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 <- ggplot(basic, aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens4 <- ggplot(basic, aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

############# polar
polar = merged_datasets %>% filter(variants %in% c("S","T","Y","N","Q"))

plot3 <- ggplot(data=polar, aes(x=met_score, y=ex14_score, color = mutation_type)) + 
  geom_point(shape=1,aes(color = mutation_type)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("polar")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 <- ggplot(polar, aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens4 <- ggplot(polar, aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

############# acidic
acidic = merged_datasets %>% filter(variants %in% c("D","E"))

plot3 <- ggplot(data=acidic, aes(x=met_score, y=ex14_score, color = mutation_type)) + 
  geom_point(shape=1,aes(color = mutation_type)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("Acidic")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 <- ggplot(acidic, aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens4 <- ggplot(acidic, aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

############# structure breaking 
structure_breaking = merged_datasets %>% filter(variants %in% c("W","Y","P","G"))

plot3 <- ggplot(data=structure_breaking, aes(x=met_score, y=ex14_score, color = mutation_type)) + 
  geom_point(shape=1,aes(color = mutation_type)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("structure_breaking")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 <- ggplot(structure_breaking, aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens4 <- ggplot(structure_breaking, aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

############# nucleophilic 
nucleophilic  = merged_datasets %>% filter(variants %in% c("S","T","C"))

plot3 <- ggplot(data=nucleophilic , aes(x=met_score, y=ex14_score, color = mutation_type)) + 
  geom_point(shape=1,aes(color = mutation_type)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("nucleophilic")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 <- ggplot(nucleophilic , aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens4 <- ggplot(nucleophilic , aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

############# amide
amide  = merged_datasets %>% filter(variants %in% c("N","Q"))

plot3 <- ggplot(data=amide , aes(x=met_score, y=ex14_score, color = mutation_type)) + 
  geom_point(shape=1,aes(color = mutation_type)) +
  theme_pubr() +
  theme(legend.position = c(0.15, 0.9))+
  xlab("MET-Ex14") + ylab("MET+Ex14")+
  ggtitle("amide")+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 <- ggplot(amide , aes(x = met_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens4 <- ggplot(nucleophilic , aes(x = ex14_score, fill = mutation_type)) + 
  geom_density(alpha = 0.4) + 
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dens3 + plot_spacer() + plot3 + dens4 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))


#####################################################################################################
### --------------------------------R and C Spine physiochemical mapping ----------------------------
#####################################################################################################

# R spine positions: L1142, M1131, F1223, H1202
# C spine positions: A1108, V1092, M1211, C1210, L1165, L2112, L1276, L1272 


#spines = scores_filtered_diff %>% filter((pos %in% c(1110, 1127,1142, 1131, 1222, 1223, 1202,1108, 1092, 1211, 1210, 1165, 2112, 1276, 1272)))
#Rspines = scores_filtered_diff %>% filter((pos %in% c(1142, 1131, 1223, 1202)))
Rspines = ex14_scores %>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1142, 1131, 1223, 1202)))
Cspines = ex14_scores %>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1108, 1092, 1211, 1212, 1210, 1165, 1276, 1272 )))


## physiochem associations 
hydrophobic<- c('I','L','V','M','A','H','K','W','Y','F')
positive<- c('R','K','H')
negative<- c('D','E')
polar_uncharged<-c('N','Q','T','S','C')

Rspine_physio<-bind_rows(
  {Rspines %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {Rspines %>% filter(variants %in% negative)%>% group_by(variants, ID='negative')},
  {Rspines %>% filter(variants %in% positive)%>% group_by(variants, ID='positive')},
  {Rspines %>% filter(variants %in% polar_uncharged)%>% group_by(variants, ID='polar uncharged')}
)

Cspine_physio<-bind_rows(
  {Cspines %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {Cspines %>% filter(variants %in% negative)%>% group_by(variants, ID='negative')},
  {Cspines %>% filter(variants %in% positive)%>% group_by(variants, ID='positive')},
  {Cspines %>% filter(variants %in% polar_uncharged)%>% group_by(variants, ID='polar uncharged')}
)

## takes the absolute sum of positional scores 
met_spines_scores <- Rspines%>% group_by(pos) %>% summarise(met_spines_abs_sum = sum(met_score, na.rm=TRUE))
ex14_Rspines_scores <- Rspines%>% group_by(pos) %>% summarise(ex14_Rspines_abs_sum= mean(ex14_score, na.rm=TRUE)) 
ex14_Cspines_scores <- Cspines%>% group_by(pos) %>% summarise(ex14_Cspines_abs_sum= mean(ex14_score, na.rm=TRUE)) 


### distributions of R and C spine as histograms 

plot_R_spine_ex14_histograms <-ggplot()+
  geom_histogram(aes(x=Rspines$IL3_withdrawal_score),color="black", fill="turquoise4", alpha = 0.5) +
  ggtitle("R-Spine Missense")+
  xlab("Activity Score")+
  ylab("Count")+
  geom_vline(xintercept=ex14_syn_mean, linetype = "dashed")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_C_spine_ex14_histograms <-ggplot()+
  geom_histogram(aes(x=Cspines$IL3_withdrawal_score),color="black", fill="seagreen4", alpha = 0.5) +
  ggtitle("C-Spine Missense")+
  xlab("Activity Score")+
  ylab("Count")+
  geom_vline(xintercept=ex14_syn_mean, linetype = "dashed")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_grid(plot_R_spine_ex14_histograms,plot_C_spine_ex14_histograms,ncol=1,nrow=2)


cbp1 <- c("#999999","#E69F00", "#56B4E9","#009E73","#E69F00",
          "#0072B2","#F0E442","#D55E00", "#CC79A7")
# violin plots with physiochem breakdown 
plot_Rspine_physio_violin <- ggplot(Rspine_physio, aes(x=factor(pos), y= IL3_withdrawal_score, fill = ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.7) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("R-spine missense substitutions")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",)
plot(plot_Rspine_physio_violin)

plot_Cspine_physio_violin <- ggplot(Cspine_physio, aes(x=factor(pos), y= IL3_withdrawal_score, fill = ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.7) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("C-spine missense substitutions")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="bottom", 
        legend.box = "horizontal",
        legend.title = element_blank())
plot(plot_Cspine_physio_violin)

plot_grid(plot_Rspine_physio_violin,plot_Cspine_physio_violin,
          nrow=2, ncol=1,rel_heights = c(1, 1.2))




################### same but for METdEx14
met_Rspines = met_scores %>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1142, 1131, 1223, 1202)))
met_Cspines = met_scores %>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1108, 1092, 1211, 1210, 1212,1165, 1276, 1272 )))


met_Rspine_physio<-bind_rows(
  {met_Rspines %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {met_Rspines %>% filter(variants %in% negative)%>% group_by(variants, ID='negative')},
  {met_Rspines %>% filter(variants %in% positive)%>% group_by(variants, ID='positive')},
  {met_Rspines %>% filter(variants %in% polar_uncharged)%>% group_by(variants, ID='polar uncharged')}
)

met_Cspine_physio<-bind_rows(
  {met_Cspines %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {met_Cspines %>% filter(variants %in% negative)%>% group_by(variants, ID='negative')},
  {met_Cspines %>% filter(variants %in% positive)%>% group_by(variants, ID='positive')},
  {met_Cspines %>% filter(variants %in% polar_uncharged)%>% group_by(variants, ID='polar uncharged')}
)

# violin plots with physiochem breakdown 
plot_met_Rspine_physio_violin <- ggplot(met_Rspine_physio, aes(x=factor(pos), y= IL3_withdrawal_score, fill = ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.6) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = met_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("R-missense substitutions")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",)
plot(plot_Rspine_physio_violin)

plot_met_Cspine_physio_violin <- ggplot(met_Cspine_physio, aes(x=factor(pos), y= IL3_withdrawal_score, fill = ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.7) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = met_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("C-spine missense substitutions")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position="bottom", 
        legend.box = "horizontal",
        legend.title = element_blank())
plot(plot_Cspine_physio_violin)

plot_grid(plot_met_Rspine_physio_violin,plot_met_Cspine_physio_violin,
          nrow=2, ncol=1,rel_heights = c(1, 1.2))



#####################################################################################################
### -------------------------------- Surface vs Core score distribution ----------------------------
#####################################################################################################

core_residues = ex14_scores%>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1095,1107,1108,1121,1124,1127,1128,1131,1139,1144,1145,
                                                                                                1146,1153,1155,1156,1161,1168,1176,1180,1181,1182,1183,
                                                                                                1184,1185,1188,1189,1191,1195,1202,1205,1206,1207,1210,
                                                                                                1212,1218,1220,1223,1226,1251,1254,1265,1266,1267,1269,
                                                                                                1270,1271,1272,1273,1275,1275,1280,1292,1312,1316,1320,
                                                                                                1321,1328,1333,1337)))

surface_residues = ex14_scores%>% filter(mutation_type != "S" & mutation_type != "N" & !(pos %in% c(1095,1107,1108,1121,1124,1127,1128,1131,1139,1144,1145,
                                                                                                1146,1153,1155,1156,1161,1168,1176,1180,1181,1182,1183,
                                                                                                1184,1185,1188,1189,1191,1195,1202,1205,1206,1207,1210,
                                                                                                1212,1218,1220,1223,1226,1251,1254,1265,1266,1267,1269,
                                                                                                1270,1271,1272,1273,1275,1275,1280,1292,1312,1316,1320,
                                                                                                1321,1328,1333,1337)))
met_core_residues = met_scores%>% filter(mutation_type != "S" & mutation_type != "N" & (pos %in% c(1095,1107,1108,1121,1124,1127,1128,1131,1139,1144,1145,
                                                                                                1146,1153,1155,1156,1161,1168,1176,1180,1181,1182,1183,
                                                                                                1184,1185,1188,1189,1191,1195,1202,1205,1206,1207,1210,
                                                                                                1212,1218,1220,1223,1226,1251,1254,1265,1266,1267,1269,
                                                                                                1270,1271,1272,1273,1275,1275,1280,1292,1312,1316,1320,
                                                                                                1321,1328,1333,1337)))

met_surface_residues = met_scores%>% filter(mutation_type != "S" & mutation_type != "N" & !(pos %in% c(1095,1107,1108,1121,1124,1127,1128,1131,1139,1144,1145,
                                                                                                    1146,1153,1155,1156,1161,1168,1176,1180,1181,1182,1183,
                                                                                                    1184,1185,1188,1189,1191,1195,1202,1205,1206,1207,1210,
                                                                                                    1212,1218,1220,1223,1226,1251,1254,1265,1266,1267,1269,
                                                                                                    1270,1271,1272,1273,1275,1275,1280,1292,1312,1316,1320,
                                                                                                    1321,1328,1333,1337)))



plot_Core_residue_hist <-ggplot()+
  geom_histogram(aes(x=core_residues$IL3_withdrawal_score),color="black", fill="indianred", alpha = 0.8) +
  ggtitle("Core missense")+
  xlab("Activity Score")+
  ylab("Count")+
  geom_vline(xintercept=ex14_syn_mean, linetype = "dashed")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_Surface_residue_hist <-ggplot()+
  geom_histogram(aes(x=surface_residues$IL3_withdrawal_score),color="black", fill="lightskyblue", alpha = 0.8) +
  ggtitle("Surface missense")+
  xlab("Activity Score")+
  ylab("Count")+
  geom_vline(xintercept=met_syn_mean, linetype = "dashed")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_grid(plot_Surface_residue_hist,plot_Core_residue_hist,ncol=2,nrow=1)


plot_MET_Core_residue_hist <-ggplot()+
  geom_histogram(aes(x=core_residues$IL3_withdrawal_score),color="black", fill="grey", alpha = 0.5) +
  geom_histogram(aes(x=met_core_residues$IL3_withdrawal_score),color="black", fill="indianred", alpha = 0.5) +
  ggtitle("Core missense")+
  xlab("Activity Score")+
  ylab("Count")+
  geom_vline(xintercept=ex14_syn_mean, linetype = "dashed")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_MET_Surface_residue_hist <-ggplot()+
  geom_histogram(aes(x=surface_residues$IL3_withdrawal_score),color="black", fill="grey", alpha = 0.5) +
  geom_histogram(aes(x=met_surface_residues$IL3_withdrawal_score),color="black", fill="lightskyblue", alpha = 0.5) +
  ggtitle("Surface missense")+
  xlab("Activity Score")+
  ylab("Count")+
  geom_vline(xintercept=met_syn_mean, linetype = "dashed")+
  theme_bw()+
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_grid(plot_MET_Surface_residue_hist,plot_MET_Core_residue_hist,ncol=1,nrow=2)


met_surface_avg_scores <- met_surface_residues %>% group_by(pos) %>% summarise(met_surface_avg = mean(IL3_withdrawal_score, na.rm=TRUE))
met_surface_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_surface_avg_3R7O, met_surface_avg_scores, "met_surface_avg")
write.pdb(x, file="met_surface_avg_3R7O.pdb")
ex14_surface_avg_scores <- surface_residues %>% group_by(pos) %>% summarise(ex14_surface_avg = mean(IL3_withdrawal_score, na.rm=TRUE))
ex14_surface_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_surface_avg_3R7O, ex14_surface_avg_scores, "ex14_surface_avg")
write.pdb(x, file="ex14_surface_avg_3R7O.pdb")



#####################################################################################################
### --------------------------------Motif Score Mapping   ----------------------------
#####################################################################################################

# JM helix : 1059, 1060, 1061, 1062, 1063, 1064,1065,1066,1067, 1068,1069,1070
# C-helix : 1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1127,1128,1129,1130,1131,1132,1133,1134
# A-loop : 1222,1223,1224,1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246
# P_loop : 1084,1085,1086,1087,1088,1089,1090,1091
# Hinge: 1158,1159,1160,1161,1162,1163,1164

JM_helix = merged_datasets %>% filter(mutation_type != "S" & (pos %in% c(1059, 1060, 1061, 1062, 1063, 1064,1065,1066,1067, 1068,1069,1070)))
C_helix = merged_datasets %>% filter(mutation_type != "S" & (pos %in% c(1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1127,1128,1129,1130,1131,1132,1133,1134 )))
A_loop = merged_datasets %>% filter(mutation_type != "S" & (pos %in% c(1222,1223,1224,1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246)))
P_loop = merged_datasets %>% filter(mutation_type != "S" & (pos %in% c(1084,1085,1086,1087,1088,1089,1090,1091)))
Hinge = merged_datasets %>% filter(mutation_type != "S" & (pos %in% c(1158,1159,1160,1161,1162,1163,1164)))


plotR1 <- ggplot(JM_helix, aes(x=factor(pos), y= met_score))+
  geom_point(shape=21,size=2)+
  #geom_boxplot()+
  ylim(-9,9)+
  xlab("Position") + ylab("Activity Score")+
  geom_hline(yintercept = -0.4341423, linetype="dashed")+
  ggtitle("METdEx14 JM-helix")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))


plotR2 <- ggplot(C_helix, aes(x=factor(pos), y= met_score))+
  geom_point(shape=21,size=2)+
  #geom_boxplot()+
  ylim(-9,9)+
  xlab("Position") + ylab("Activity Score")+
  geom_hline(yintercept = -0.4341423, linetype="dashed")+
  ggtitle("METdEx14 C_helix")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))


plotR3 <- ggplot(A_loop, aes(x=factor(pos), y= met_score))+
  geom_point(shape=21,size=2)+
  #geom_boxplot()+
  ylim(-9,9)+
  xlab("Position") + ylab("Activity Score")+
  geom_hline(yintercept = -0.4341423, linetype="dashed")+
  ggtitle("METdEx14 A_loop")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

plotR4 <- ggplot(P_loop, aes(x=factor(pos), y= met_score))+
  geom_point(shape=21,size=2)+
  #geom_boxplot()+
  ylim(-9,9)+
  xlab("Position") + ylab("Activity Score")+
  geom_hline(yintercept = -0.4341423, linetype="dashed")+
  ggtitle("METdEx14 P_loop")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

plotR5 <- ggplot(Hinge, aes(x=factor(pos), y= met_score))+
  geom_point(shape=21,size=2)+
  #geom_boxplot()+
  ylim(-9,9)+
  xlab("Position") + ylab("Activity Score")+
  geom_hline(yintercept = -0.4341423, linetype="dashed")+
  ggtitle("METdEx14 Hinge")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

plotR6 <- ggplot(JM_helix, aes(x=ex14_score))+
  geom_histogram(color="black",fill="darkseagreen3", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = -0.1426368, linetype="dashed")+
  ggtitle("JM-helix")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

plotR7 <- ggplot(C_helix, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = -0.1426368, linetype="dashed")+
  ggtitle("C-helix")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot_grid(plotR6, plotR7, ncol=2, nrow=1)


plotR8 <- ggplot(A_loop, aes(x=factor(pos), y= ex14_score))+
  geom_point(shape=21,size=2)+
  #geom_boxplot()+
  ylim(-9,9)+
  xlab("Position") + ylab("Activity Score")+
  geom_hline(yintercept = -0.1426368, linetype="dashed")+
  ggtitle("MET+Ex14 A_loop")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
plotR9 <- ggplot(P_loop, aes(x=factor(pos), y= ex14_score))+
  geom_point(shape=21,size=2)+
  #geom_boxplot()+
  ylim(-9,9)+
  xlab("Position") + ylab("Activity Score")+
  geom_hline(yintercept = -0.1426368, linetype="dashed")+
  ggtitle("MET+Ex14 P_loop")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

plotR10 <- ggplot(Hinge, aes(x=ex14_score))+
  geom_histogram()+
  #geom_boxplot()+
  xlab("Activity Score") + ylab("Frequency")+
  xlim(-10,10)+
  geom_vline(xintercept = -0.1426368, linetype="dashed")+
  ggtitle("MET+Ex14 Hinge")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
plot(plotR10)

P1153 = merged_datasets %>% filter(mutation_type != "N" & mutation_type != "S" & (pos %in% c(1153)))

plot_P1153 <- ggplot()+
  geom_histogram(aes(x=P1153$ex14_score),color="black",fill="grey", alpha=0.5)+ 
  geom_histogram(aes(x=P1153$met_score),color="black",fill="deepskyblue2", alpha=0.5)+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed", color="grey")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed", color="deepskyblue2")+
  ggtitle("P1153")+
  theme_bw() +
  theme(
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_P1153 )


plot_grid(plotR1,plotR2,plotR3,plotR4,plotR5,
          plotR6,plotR7,plotR8,plotR9,plotR10,
          ncol=5,nrow=2)


####################################################################
####----------------------- interface analysis ---------------------
####################################################################

# here there are some important residues within the parallel interface 
# I want to dissect any interactions there that could support a differential 
# the positions are H1068, E1061 which display a salt bridge in the parallel interface 
# I also want to look at I1130, S1122,I1118

JM_interface = merged_datasets %>% filter(mutation_type != "N" & mutation_type != "S" &  (pos %in% c(1061,1065,1068)))
Chelix_interface =  merged_datasets %>% filter(mutation_type != "N" & mutation_type != "S" & (pos %in% c(1118,1122,1126,1130)))

#Chelix_surface = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1118,1122,1130)))


H1068 = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1068)))
E1061 = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1061)))

I1130 = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1130)))
S1122 = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1122)))
I1118 = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1118)))
T1126 = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1126)))
Q1123 = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1123)))



############ violin plots of surface residues of the MET+/-Ex14 data

hydrophobic<- c('I','L','V','M','A','W','Y','F')
positive<- c('R','K','H')
negative<- c('D','E')
polar_uncharged<-c('N','Q','T','S','C')

JM_surface_physio<-bind_rows(
  {JM_interface %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {JM_interface %>% filter(variants %in% negative)%>% group_by(variants, ID='negative')},
  {JM_interface %>% filter(variants %in% positive)%>% group_by(variants, ID='positive')},
  {JM_interface %>% filter(variants %in% polar_uncharged)%>% group_by(variants, ID='polar uncharged')}
)

Chelix_surface_physio<-bind_rows(
  {Chelix_interface %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {Chelix_interface %>% filter(variants %in% negative)%>% group_by(variants, ID='negative')},
  {Chelix_interface %>% filter(variants %in% positive)%>% group_by(variants, ID='positive')},
  {Chelix_interface %>% filter(variants %in% polar_uncharged)%>% group_by(variants, ID='polar uncharged')}
)



ex14_JM_interface <- ggplot(JM_surface_physio, aes(x=factor(pos), y= ex14_score, fill =ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.5) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle(expression("MET+Ex14 JM-helix"))+
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none")
plot(ex14_JM_interface)

ex14_Chelix_interface <- ggplot(Chelix_surface_physio, aes(x=factor(pos), y= ex14_score, fill = ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.5) + 
  ylim(-8,1)+
  scale_fill_manual(values=cbp1)+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle(expression("MET+Ex14 C-helix"))+
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none")
plot(ex14_Chelix_interface)


met_JM_interface <- ggplot(JM_surface_physio, aes(x=factor(pos), y= met_score, fill = ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.5) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = met_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("METdEx14 JM-helix")+
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none")
plot(met_JM_interface)

met_Chelix_interface <- ggplot(Chelix_surface_physio, aes(x=factor(pos), y= met_score, fill = ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.5) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = met_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle("METdEx14 C-helix")+
  theme_bw()+
  theme(text = element_text(size = 12),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(met_Chelix_interface)

plot_grid(ex14_JM_interface,met_JM_interface,
          ex14_Chelix_interface,met_Chelix_interface,
          ncol=4,nrow=,rel_widths = c(1,1,1,1.7))



###################### interface comparisons

parallel_interface = merged_datasets %>% filter(mutation_type != "N" & mutation_type != "S" & (pos %in% c(1061,1065,1068,1118,1122,1130)))
head_to_tail_interface = merged_datasets %>% filter(mutation_type != "N" & mutation_type != "S"  & (pos %in% c(1061,1064,1068,1067,1340,1336,1343,1314)))
head_to_head_interface = merged_datasets %>% filter(mutation_type != "N" & mutation_type != "S" & (pos %in% c(1060,1063,1064,1067,1070,1114,1113,1116)))

plot_parallel_interface <- ggplot()+
  geom_point(aes(x=head_to_tail_interface$ex14_score,y=head_to_tail_interface$met_score,color="Head-to-tail"),show.legend = T)+
  geom_point(aes(x=head_to_head_interface$ex14_score,y=head_to_head_interface$met_score,color="Head-to-head"),show.legend = T)+
  geom_point(aes(x=parallel_interface$ex14_score,y=parallel_interface$met_score,color="Parallel"),show.legend = T)+
  geom_vline(xintercept = -2.5,linetype="dashed")+
  geom_hline(yintercept = -2.5,linetype="dashed")+
  ylim(-10,3)+
  xlim(-10,3)+
  theme_bw()+
  theme(text = element_text(size = 12),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())+
  xlab("MET+Ex14 Activity Score") + ylab("METdEx14 Activity Score")+
  scale_colour_manual(name = 'Interface', 
                      guide = 'legend',
                      values = c('Parallel'='black','Head-to-tail'='#E69F00',
                                 'Head-to-head'='#56B4E9'), 
                      labels = c('Parallel',
                                 'Head-to-tail',
                                 'Head-to-head '))
plot(plot_parallel_interface)

###################### JM surface vs core comparisions 


JM_surface = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(11059,1060,1061,1064,1065,1067,1068,1070)))
JM_core = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1062,1063,1066,1067,1069)))
Chelix_surface = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1117,1120,1123,1126,1132,1133,1334)))
Chelix_core = merged_datasets %>% filter(mutation_type != "N" & (pos %in% c(1118,1119,1121,1122,1124,1125,1127,1128,1129,1131)))


ex14_JMsurface <- ggplot(JM_surface, aes(x=ex14_score))+
  geom_histogram(color="black",fill="orange", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle(expression("MET+Ex14 JM-helix surface"))+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

ex14_Chelix_surface <- ggplot(Chelix_surface, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 C-helix surface")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

met_JMsurface  <- ggplot(JM_surface , aes(x=met_score))+
  geom_histogram(color="black",fill="orange", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 JM-helix surface")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

met_Chelix_surface  <- ggplot(Chelix_surface , aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 C-helix surfaces")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())


plot_grid(ex14_JMsurface,ex14_Chelix_surface,
          met_JMsurface,met_Chelix_surface, 
          nrow=2,ncol=2)

plot_JM_surface_comparisons <-ggplot(JM_surface)+
   geom_histogram(aes(x=ex14_score, color = "#F39C12", fill="MET+Ex14"), alpha = 0.8) + 
   geom_histogram(aes(x=met_score, color = "#ABB2B9", fill="METdEx14"), alpha = 0.3) +
   scale_fill_manual(name="",values = c("#ABB2B9","#F39C12"))+
   scale_colour_manual(name="",values = c("#F39C12","#ABB2B9"), guide=FALSE)+
   ggtitle("JM-helix surface")+
   xlab("Activity Score")+
   ylab("Count")+
   theme_bw() +
   theme(text = element_text(size = 15),
         panel.background=element_rect(colour="black",size=1),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.background = element_blank(),
         legend.position = c(0.3,.75))

plot_JM_core_comparisons <-ggplot(JM_core)+
  geom_histogram(aes(x=ex14_score, color = "#F39C12", fill="MET+Ex14"), alpha = 0.8) + 
  geom_histogram(aes(x=met_score, color = "#ABB2B9", fill="METdEx14"), alpha = 0.3) +
  scale_fill_manual(name="",values = c("#ABB2B9","#F39C12"))+
  scale_colour_manual(name="",values = c("#F39C12","#ABB2B9"), guide=FALSE)+
  ggtitle("JM-helix core")+
  xlab("Activity Score")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = c(0.3,.75))

 plot_grid(plot_JM_surface_comparisons,plot_JM_core_comparisons,
           ncol=2,nrow=1,
           rel_widths = c(1,1))
 
 plot_Chelix_surface_comparisons <-ggplot(Chelix_surface)+
   geom_histogram(aes(x=ex14_score, color = "#F39C12", fill="MET+Ex14"), alpha = 0.8) + 
   geom_histogram(aes(x=met_score, color = "#ABB2B9", fill="METdEx14"), alpha = 0.3) +
   scale_fill_manual(name="",values = c("#ABB2B9","#F39C12"))+
   scale_colour_manual(name="",values = c("#F39C12","#ABB2B9"), guide=FALSE)+
   ggtitle("C-helix surface")+
   xlab("Activity Score")+
   ylab("Count")+
   theme_bw() +
   theme(text = element_text(size = 15),
         panel.background=element_rect(colour="black",size=1),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.background = element_blank(),
         legend.position = c(0.2,.85))
 
 plot_Chelix_core_comparisons <-ggplot(Chelix_core)+
   geom_histogram(aes(x=ex14_score, color = "#F39C12", fill="MET+Ex14"), alpha = 0.8) + 
   geom_histogram(aes(x=met_score, color = "#ABB2B9", fill="METdEx14"), alpha = 0.3) +
   scale_fill_manual(name="",values = c("#ABB2B9","#F39C12"))+
   scale_colour_manual(name="",values = c("#F39C12","#ABB2B9"), guide=FALSE)+
   ggtitle("C-helix core")+
   xlab("Activity Score")+
   theme_bw() +
   theme(text = element_text(size = 15),
         panel.background=element_rect(colour="black",size=1),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.background = element_blank(),
         legend.position = c(0.2,.85))
 
 plot_grid(plot_Chelix_surface_comparisons,plot_Chelix_core_comparisons,
           ncol=2,nrow=1,
           rel_widths = c(1,1))
 




#### core vs surface 

ex14_Chelix_core <- ggplot(Chelix_core, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 C-helix core")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

met_Chelix_core  <- ggplot(Chelix_core , aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 C-helix core")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

plot_grid(ex14_Chelix_core,ex14_Chelix_surface,
          met_Chelix_core,met_Chelix_surface, 
          nrow=2,ncol=2)



ex14_JM_core <- ggplot(JM_core, aes(x=ex14_score))+
  geom_histogram(color="black",fill="orange", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 JM-helix core")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

met_JM_core  <- ggplot(JM_core , aes(x=met_score))+
  geom_histogram(color="black",fill="orange", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 JM-helix core")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

plot_grid(ex14_JM_core,ex14_JMsurface,
          met_JM_core,met_JMsurface, 
          ex14_Chelix_core,ex14_Chelix_surface,
          met_Chelix_core,met_Chelix_surface, 
          nrow=4,ncol=2)







#### specific residues

ex14_H1068 <- ggplot(H1068, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 H1068")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(ex14_H1068)


ex14_E1061 <- ggplot(E1061, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 E1061")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(ex14_E1061 )

met_H1068  <- ggplot(H1068, aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 H1068")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(met_H1068)

met_E1061  <- ggplot(E1061, aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 E1061")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(met_E1061)

plot_grid(ex14_H1068,ex14_E1061,met_H1068,met_E1061,nrow=2,ncol=2)


ex14_I1130 <- ggplot(I1130, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 I1130")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

ex14_I1118 <- ggplot(I1118, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 I1118")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

ex14_S1122 <- ggplot(S1122, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 S1122")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

ex14_T1126 <- ggplot(T1126, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 T1126")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

ex14_Q1123 <- ggplot(Q1123, aes(x=ex14_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = ex14_syn_mean, linetype="dashed")+
  ggtitle("MET+Ex14 Q1123")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())




met_I1130  <- ggplot(I1130 , aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 I1130")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

met_I1118  <- ggplot(I1118 , aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 I1118")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

met_S1122  <- ggplot(S1122, aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 S1122")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

met_T1126  <- ggplot(T1126, aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 T1126")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())

met_Q1123  <- ggplot(Q1123, aes(x=met_score))+
  geom_histogram(color="black",fill="deepskyblue2", alpha=0.5)+
  #geom_boxplot()+
  xlim(-10,10)+
  xlab("Activity Score") + ylab("Count")+
  geom_vline(xintercept = met_syn_mean, linetype="dashed")+
  ggtitle("METdEx14 Q1123")+
  theme_bw() +
  theme(text = element_text(size = 15),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())


plot_grid(ex14_I1130,ex14_I1118,ex14_S1122,ex14_T1126,ex14_Q1123,
          met_I1130,met_I1118,met_S1122,met_T1126,met_Q1123,
          nrow=2,ncol=5)


################################################################################################
#------------------JM helix Hydrophobic Staple  positional analysis- (Fig 3)--------------------
################################################################################################
## try to fill with gradient of color based on physiochem
## hydrophobic, negative, positive, aromatic, small, helix break 


JM_helix = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N"&(pos %in% c(1062,1063,1066,1069,1070)))
C_helix = merged_datasets %>% filter(mutation_type != "S" & mutation_type != "N" &(pos %in% c(1118,1121,1122,1125,1129)))


hydrophobic<- c('I','L','V','M','A','W','Y','F','C')
not_hyrdrophobic<-c('N','Q','T','S','D','E','R','K','H')
helix_break <- c('G','P')

JM_physio<-bind_rows(
  {JM_helix %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {JM_helix%>% filter(variants %in% not_hyrdrophobic)%>% group_by(variants, ID='not_hyrdrophobic')},
  {JM_helix%>% filter(variants %in% helix_break)%>% group_by(variants, ID='helix_break')}
)

Chelix_physio<-bind_rows(
  {C_helix %>% filter(variants %in% hydrophobic)%>% group_by(variants, ID='hydrophobic')},
  {C_helix %>% filter(variants %in% not_hyrdrophobic)%>% group_by(variants, ID='not_hyrdrophobic')},
  {C_helix%>% filter(variants %in% helix_break)%>% group_by(variants, ID='helix_break')}
)


plot_JM_helix <- ggplot(JM_physio, aes(x=factor(pos), y= ex14_score, fill =ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.3) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle(expression("JM-helix missense "))+
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "none")
plot(plot_JM_helix)

plot_C_helix <- ggplot(Chelix_physio, aes(x=factor(pos), y= ex14_score, fill =ID))+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 1.3) + 
  scale_fill_manual(values=cbp1)+
  ylim(-8,1)+
  geom_hline(yintercept = ex14_syn_mean, linetype="dashed")+
  xlab("Position") + ylab("Activity Score")+
  ggtitle(expression("C-helix missense"))+
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.background=element_rect(colour="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank())
plot(plot_C_helix)

plot_grid(plot_JM_helix,plot_C_helix,ncol=2,nrow=1, rel_widths = c(1,1.5))


### N-lobe hydrophobic network 
N_hydro_network = merged_datasets %>% filter(mutation_type != "S" & (pos %in% c(1071,1076,1097,1156)))
plot<-ggplot(N_hydro_network, aes(x=factor(pos), y= ex14_score, fill = factor(pos)))+
  geom_violin(alpha = 0.2) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5) + 
  scale_fill_brewer(palette = "Oranges")+
  ylim(-10,5)+
  geom_hline(yintercept = -0.1426368, linetype="dashed")+
  ggtitle("N-lobe Hydrophobic Network")+
  xlab("Position")+
  ylab("Activty Score")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text=element_text(size=15),
        legend.position = "none")
plot(plot)


################################################################################################
#------------------ RMSD analysis of active and inactive structures (fig 3) -------------------
################################################################################################
# utilizes the bio3D package for RMSD analysis

MET_active<- read.pdb("3R7O")
MET_active_dcd <- read.dcd("3R7O")

MET_inactive<- read.pdb("5HTI")

ca.met_active.inds <- atom.select(MET_active, elety="CA")
MET_active_xyz<- fit.xyz(fixed=MET_active$xyz, mobile=MET_active_dcd,
                           fixed.inds=ca.met_active.inds$xyz,
                           mobile.inds=ca.met_active.inds$xyz)

rd <- rmsd(MET_active_xyz[1,ca.met_active.inds$xyz], ca.met_active.inds[,ca.met_active.inds$xyz])

plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)



###--------------------- plot of hydrophobic scores----------
met_spines_hydrophobic = met_spines %>% filter(variants %in% c("V","L","I","M","P","G","A","W","F"))
ex14_spines_hydrophobic = ex14_spines %>% filter(variants %in% c("V","L","I","M","P","G","A","W","F"))

met_spines_hydrophobic_avg_scores <- met_spines_hydrophobic %>% group_by(pos) %>% summarise(met_spines_hydrophobic_avg = mean(IL3_withdrawal_score, na.rm=TRUE))
ex14_spines_hydrophobic_avg_scores <- ex14_spines_hydrophobic %>% group_by(pos) %>% summarise(ex14_spines_hydrophobic_avg = mean(IL3_withdrawal_score, na.rm=TRUE))   

met_spine_hydrophobic_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_spine_hydrophobic_avg_3R7O, met_spines_hydrophobic_avg_scores, "met_spines_hydrophobic_avg")
write.pdb(x, file="met_spine_hydrophobic_avg_3R7O.pdb")
ex14_spines_hydrophobic_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_spines_hydrophobic_avg_3R7O, ex14_spines_hydrophobic_avg_scores , "ex14_spines_hydrophobic_avg")
write.pdb(x, file="ex14_spine_hydrophobic_avg_3R7O.pdb")



###----------------------- plot of aromatic scores
met_spines_aromatic = met_spines %>% filter(variants %in% c("Y","W","F"))
ex14_spines_aromatic = ex14_spines %>% filter(variants %in% c("Y","W","F"))

met_spines_aromatic_avg_scores <- met_spines_aromatic %>% group_by(pos) %>% summarise(met_spines_aromatic_avg = mean(IL3_withdrawal_score, na.rm=TRUE))
ex14_spines_aromatic_avg_scores <- ex14_spines_aromatic %>% group_by(pos) %>% summarise(ex14_spines_aromatic_avg = mean(IL3_withdrawal_score, na.rm=TRUE)) 

met_spine_aromatic_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_spine_aromatic_avg_3R7O, met_spines_aromatic_avg_scores, "met_spines_aromatic_avg")
write.pdb(x, file="met_spine_aromatic_avg_3R7O.pdb")

ex14_spines_aromatic_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_spines_aromatic_avg_3R7O, ex14_spines_aromatic_avg_scores, "ex14_spines_aromatic_avg")
write.pdb(x, file="ex14_spine_aromatic_avg_3R7O.pdb")


###----------------------- plot of acidic scores
met_spines_acidic = met_spines %>% filter(variants %in% c("D","E"))
ex14_spines_acidic = ex14_spines %>% filter(variants %in% c("D","E"))

met_spines_acidic_avg_scores <- met_spines_acidic %>% group_by(pos) %>% summarise(met_spines_acidic_avg = mean(IL3_withdrawal_score, na.rm=TRUE))
ex14_spines_acidic_avg_scores <- ex14_spines_acidic %>% group_by(pos) %>% summarise(ex14_spines_acidic_avg = mean(IL3_withdrawal_score, na.rm=TRUE)) 

met_spine_acidic_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_spine_acidic_avg_3R7O, met_spines_acidic_avg_scores, "met_spines_acidic_avg")
write.pdb(x, file="met_spine_acidic_avg_3R7O.pdb")

ex14_spines_acidic_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(ex14_spines_acidic_avg_3R7O, ex14_spines_acidic_avg_scores , "ex14_spines_acidic_avg")
write.pdb(x, file="ex14_spine_acidic_avg_3R7O.pdb")


###----------------------- plot basic scores 
met_spines_basic = met_spines %>% filter(variants %in% c("H","K","R"))
ex14_spines_basic = ex14_spines %>% filter(variants %in% c("H","K","R"))

met_spines_basic_avg_scores <- met_spines_basic %>% group_by(pos) %>% summarise(met_spines_basic_avg = mean(IL3_withdrawal_score, na.rm=TRUE))
ex14_spines_basic_avg_scores <- ex14_spines_basic %>% group_by(pos) %>% summarise(ex14_spines_basic_avg = mean(IL3_withdrawal_score, na.rm=TRUE)) 

met_spine_basic_avg_3R7O <- read.pdb("3R7O")
x = map_scores_pdb(met_spine_basic_avg_3R7O, met_spines_basic_avg_scores, "met_spines_basic_avg")
write.pdb(x, file="met_spine_basic_avg_3R7O.pdb")


############################# Protein Structure Mapping #####################################################

###----------- map GOF and LOF mutations on PDBs  ----------------------

### color with:  spectrum b, white_sky_blue, minimum = 0, maximum=5
### average GOF on active structures 
met_GOF_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_GOF_avg_3R7O, GOF_avg_df, "GOF_met")
write.pdb(x, file="met_GOF_avg_3R7O.pdb")

ex14_GOF_avg_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_GOF_avg_3R7O, GOF_avg_df, "GOF_ex14")
write.pdb(x, file="ex14_GOF_avg_3R7O.pdb")

### max GOF on active structures 
met_GOF_max_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_GOF_max_3R7O, GOF_max_df, "max_GOF_met")
write.pdb(x, file="met_GOF_max_3R7O.pdb")

ex14_GOF_max_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_GOF_max_3R7O, GOF_max_df, "max_GOF_ex14")
write.pdb(x, file="ex14_GOF_max_3R7O.pdb")

### max GOF on inactive structures 
met_GOF_max_2G15 <- read.pdb("2G15")
x = map_scores_pdb(met_GOF_max_2G15, GOF_max_df, "max_GOF_met")
write.pdb(x, file="met_GOF_max_2G15.pdb")

ex14_GOF_max_2G15<- read.pdb("2G15")
x = map_scores_pdb(ex14_GOF_max_2G15, GOF_max_df, "max_GOF_ex14")
write.pdb(x, file="ex14_GOF_max_2G15.pdb")


### color with: spectrum b, firebrick_white, minimum = -15, maximum=-5

### avg LOF on inactive structures 
met_LOF_avg_2G15<- read.pdb("2G15")
x = map_scores_pdb(met_LOF_avg_2G15, LOF_avg_df, "LOF_met")
write.pdb(x, file="met_LOF_avg_2G15.pdb")

ex14_LOF_avg_2G15<- read.pdb("2G15")
x = map_scores_pdb(ex14_LOF_avg_2G15, LOF_avg_df, "LOF_ex14")
write.pdb(x, file="ex14_LOF_avg_2G15.pdb")


### min LOF on inactive structures 
met_LOF_min_2G15<- read.pdb("2G15")
x = map_scores_pdb(met_LOF_min_2G15, LOF_min_df, "min_LOF_met")
write.pdb(x, file="met_LOF_min_2G15.pdb")

ex14_LOF_min_2G15<- read.pdb("2G15")
x = map_scores_pdb(ex14_LOF_min_2G15, LOF_min_df, "min_LOF_ex14")
write.pdb(x, file="ex14_LOF_min_2G15.pdb")

### min LOF on inactive structures 
met_LOF_min_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(met_LOF_min_3R7O, LOF_min_df, "min_LOF_met")
write.pdb(x, file="met_LOF_min_3R7O.pdb")

ex14_LOF_min_3R7O<- read.pdb("3R7O")
x = map_scores_pdb(ex14_LOF_min_3R7O, LOF_min_df, "min_LOF_ex14")
write.pdb(x, file="ex14_LOF_min_3R7O.pdb")




###########################

### Physicochemistry heatmaps

# hydrophobic: V,L,I,M,P,G,A,W,F
# aromatic : F,Y,W 
# acidic : D,E 
# basic : H,K,R
# amide : N,Q 
# nucleophilic: S,T,C
# polar, uncharged: S,T,Y,N,Q, C
# small: G, A
# structure breaking: P,G
# aliphatic: I,L,V,A
# phosphorylation : S,T,Y


order <- c('aromatic','aliphatic', 'hydrophobic','branched', 'tiny', 'small',
           'polar','charged','negative','positive','proton_don_acc','helix_break','phosphorylation_comp')
hydrophobic<- c('I','L','V','M','C','A','T','H','K','W','Y','F')
aromatic<- c('W','F','Y','H')
aliphatic<- c('V','A','L','I')
small<- c('C','V','A','T','S','G')
tiny<- c('G','S','C','A')
positive<- c('R','K','H')
negative<- c('D','E')
charged<- c('R','K','H','D','E')
polar<- c('N','Q','C','S','H','Y','T','R','K','D','E')
branched<-c('V','I','T')
helix_break<- c('G','P')
proton_don_acc<-c('C','H')
phosphorylation_comp <- c('S','T','Y')

met_wt_1 = str_split(substr(met_wt_sequence, 1, 100), '')[[1]]
met_wt_2 = str_split(substr(met_wt_sequence, 101, 200), '')[[1]]
met_wt_3 = str_split(substr(met_wt_sequence, 201, 287), '')[[1]]


physicochemistry_means<-bind_rows(
  {met_scores_ %>% filter(variants %in% helix_break)%>% group_by(pos) %>% summarize(ID='helix_break',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% hydrophobic)%>% group_by(pos) %>% summarize(ID='hydrophobic',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% aliphatic)%>% group_by(pos) %>% summarize(ID='aliphatic',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% aromatic)%>% group_by(pos) %>% summarize(ID='aromatic',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% negative)%>% group_by(pos) %>% summarize(ID='negative',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% positive)%>% group_by(pos) %>% summarize(ID='positive',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% tiny)%>% group_by(pos) %>% summarize(ID='tiny',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% small)%>% group_by(pos) %>% summarize(ID='small',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% charged)%>% group_by(pos) %>% summarize(ID='charged',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% polar)%>% group_by(pos) %>% summarize(ID='polar',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% phosphorylation_comp)%>% group_by(pos) %>% summarize(ID='phosphorylation_comp',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% branched)%>% group_by(pos) %>% summarize(ID='branched',MET_score=mean(IL3_withdrawal_score))},
  {met_scores_ %>% filter(variants %in% proton_don_acc)%>% group_by(pos) %>% summarize(ID='proton_don_acc',MET_score=mean(IL3_withdrawal_score))},
)

met_row4 =  ggplot(data = physicochemistry_means[physicochemistry_means$pos %in% c(1059:1158),], 
         aes(x = pos, y = factor(ID, level = order), fill = MET_score)) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,3)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1158, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1158),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")


met_row5 = ggplot(data = physicochemistry_means[physicochemistry_means$pos %in% c(1159:1258),], 
        aes(x = pos, y = factor(ID, level = order), fill = MET_score)) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,3)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1159,1258, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1159,1258),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")


met_row6 = ggplot(data = physicochemistry_means[physicochemistry_means$pos %in% c(1259:1345),], 
        aes(x = pos, y = factor(ID, level = order), fill = MET_score)) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,3)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1259,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1259,1345),
                       labels = met_wt_3,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")


met_DMS = ggarrange(met_row4, met_row5, met_row6,
                nrow =3, ncol = 1)

ggsave("Plots/met_physiochem_heatmap.pdf", height = 7, width = 8.5, met_DMS)
ggsave("Plots/met_physiochem_heatmap.png", height = 7, width = 8.5, met_DMS)


########################


physicochemistry_means_ex14<-bind_rows(
    {ex14_scores_ %>% filter(variants %in% helix_break)%>% group_by(pos) %>% summarize(ID='helix_break',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% hydrophobic)%>% group_by(pos) %>% summarize(ID='hydrophobic',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% aliphatic)%>% group_by(pos) %>% summarize(ID='aliphatic',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% aromatic)%>% group_by(pos) %>% summarize(ID='aromatic',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_%>% filter(variants %in% negative)%>% group_by(pos) %>% summarize(ID='negative',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% positive)%>% group_by(pos) %>% summarize(ID='positive',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% tiny)%>% group_by(pos) %>% summarize(ID='tiny',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% small)%>% group_by(pos) %>% summarize(ID='small',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% charged)%>% group_by(pos) %>% summarize(ID='charged',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% polar)%>% group_by(pos) %>% summarize(ID='polar',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% phosphorylation_comp)%>% group_by(pos) %>% summarize(ID='phosphorylation_comp',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% branched)%>% group_by(pos) %>% summarize(ID='branched',Ex14_score=mean(IL3_withdrawal_score))},
    {ex14_scores_ %>% filter(variants %in% proton_don_acc)%>% group_by(pos) %>% summarize(ID='proton_don_acc',Ex14_score=mean(IL3_withdrawal_score))},
  )

ex14_row4 =  ggplot(data = physicochemistry_means_ex14[physicochemistry_means_ex14$pos %in% c(1059:1158),], 
                   aes(x = pos, y = factor(ID, level = order), fill = Ex14_score)) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,2)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1059,1158, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1059,1158),
                       labels = met_wt_1,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")


ex14_row5 = ggplot(data = physicochemistry_means_ex14[physicochemistry_means_ex14$pos %in% c(1159:1258),], 
                  aes(x = pos, y = factor(ID, level = order), fill = Ex14_score)) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,2)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1159,1258, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1159,1258),
                       labels = met_wt_2,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

ex14_row6 = ggplot(data = physicochemistry_means_ex14[physicochemistry_means_ex14$pos %in% c(1259:1345),], 
                  aes(x = pos, y = factor(ID, level = order), fill = Ex14_score)) +
  geom_tile(size = 0.2) +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,l1=0.2, l3 = 0.2, p1=0.9, p3 = .4, p4 = 0.9,rev=FALSE, na.value = 'lightyellow',limits=c(-10,2)) + 
  scale_color_manual(values = c(NA,'green')) +
  scale_x_continuous(breaks = seq(1259,1345, by = 5),
                     expand = c(0,0),
                     sec.axis = sec_axis(
                       trans = ~.,
                       name = "Sequence",
                       breaks = seq(1259,1345),
                       labels = met_wt_3,
                       guide = derive()
                     )) +
  coord_fixed(ratio = 1) +
  theme(
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "lightgray"), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.text = element_text(size = 4),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, family = "mono"),
    axis.ticks = element_blank(),
    legend.position="none"
  ) +
  labs(y = "Mutation", x = "Position")

met_DMS = ggarrange(ex14_row4, ex14_row5, ex14_row6,
                    nrow =3, ncol = 1)

ggsave("Plots/ex14_physiochem_heatmap.pdf", height = 7, width = 8.5, met_DMS)
ggsave("Plots/ex14_physiochem_heatmap.png", height = 7, width = 8.5, met_DMS)






