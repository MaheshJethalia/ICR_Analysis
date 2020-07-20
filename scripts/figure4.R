library(data.table)
library(ggplot2)
library(doMC)
library(heatmap3)
library(gplots)
library(igraph)
library(lattice)
library(parallel)
library(reshape2)
#source("https://bioconductor.org/biocLite.R")
#biocLite('org.Hs.eg.db')
library('org.Hs.eg.db')
#biocLite("fgsea")
require(fgsea)
library(fastcluster)
library(factoextra)
library(devtools)
library(imager)
library(gridExtra)
library(purrr)
library(dplyr)
library(jpeg)
library(png)
library(grid)
library(plyr)
library(RColorBrewer)
library(ggsignif)
library(plotly)
library(processx)
library(ggpubr)
library(grImport2)
library(clusterProfiler)
library(magrittr)
library(msigdbr)
library(ComplexHeatmap)
library(networkD3)
warnings("off")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
#source("scripts/heatmap.3.R")

registerDoMC(20)

setwd("../ICR_All_Info/")

source('scripts/mra-analysis.R')
source('scripts/gene-reverse-network.R')
source('scripts/get_functions.R')

icr_enabled <- c("BLCA","BRCA","HNSC","LIHC","SARC","SKCM","STAD","UCEC")
icr_disabled <- c("LGG","KIRC","PAAD","UVM")

common_tfs_enabled_list <- NULL
common_tfs_disabled_list <- NULL

output_df <- read.table("Results/Text_Results/All_TF_Activity_Information.csv",header = TRUE)
output_df$TF <- as.character(as.vector(output_df$TF))
output_df$Cancer <- as.character(as.vector(output_df$Cancer))
colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#a9a9a9', '#008080', '#e6beff')

common_tfs_enabled_list <- c(output_df[output_df$Cancer==icr_enabled[1],]$TF)
for (i in 2:length(icr_enabled))
{
  new_tfs <- output_df[output_df$Cancer==icr_enabled[i],]$TF
  common_tfs_enabled_list <- intersect(common_tfs_enabled_list,new_tfs)
}

first_cancer <- icr_enabled[1];
load(paste0("Results/",first_cancer,"/Adjacency_Matrix/",first_cancer,"_Full_TopMR_Info_FGSEA_BC.Rdata"))
topmr <- topmr_info$pathway
list_common_mrs_enabled <- as.character(as.vector(topmr));
topMRs_enabled <- rep(list("BLCA"),length(list_common_mrs_enabled))
names(topMRs_enabled) <- list_common_mrs_enabled
for(i in 2:length(icr_enabled))
{
  cancer_type <- icr_enabled[i]
  load(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_TopMR_Info_FGSEA_BC.Rdata"))
  topmr <- topmr_info$pathway
  list_common_mrs_enabled <- as.character(as.vector(topmr));
  for (j in 1:length(list_common_mrs_enabled))
  {
    if (list_common_mrs_enabled[j] %in% names(topMRs_enabled))
    {
      topMRs_enabled[[list_common_mrs_enabled[j]]] <- c(topMRs_enabled[[list_common_mrs_enabled[j]]],cancer_type)
    } else {
      topMRs_enabled[[list_common_mrs_enabled[j]]] <- cancer_type
    }
  }
}

all_icr_cancers <- unlist(lapply(topMRs_enabled,function(x) paste(x,collapse=" ")))
no_icr_cancers <- unlist(lapply(topMRs_enabled,function(x) length(x)))
enabled_df <- cbind(names(all_icr_cancers),as.character(all_icr_cancers),as.numeric(no_icr_cancers))
enabled_df <- as.data.frame(enabled_df)
colnames(enabled_df) <- c("MR","List_ICR_Cancers","No_ICR_Cancers")
enabled_df$MR <- as.character(as.vector(enabled_df$MR))
enabled_df$List_ICR_Cancers <- as.character(as.vector(enabled_df$List_ICR_Cancers))
enabled_df$No_ICR_Cancers <- as.character(as.vector(enabled_df$No_ICR_Cancers))
enabled_df <- enabled_df[order(enabled_df$No_ICR_Cancers,decreasing=T),]

enabled_8_MR_and_TF <- enabled_df[enabled_df$No_ICR_Cancers==8,]$MR
enabled_7_MR <- enabled_df[enabled_df$No_ICR_Cancers==7,]$MR
enabled_7_MR_and_TF <- enabled_7_MR[enabled_7_MR %in% common_tfs_enabled_list]
enabled_7_MR_and_TF <- paste0(enabled_7_MR_and_TF,collapse = ", ")
enabled_6_MR <- enabled_df[enabled_df$No_ICR_Cancers==6,]$MR
enabled_6_MR_and_TF <- enabled_6_MR[enabled_6_MR %in% common_tfs_enabled_list]
enabled_6_MR_and_TF <- paste0(enabled_6_MR_and_TF,collapse = ", ")
enabled_5_MR <- enabled_df[enabled_df$No_ICR_Cancers==5,]$MR
enabled_5_MR_and_TF <- enabled_5_MR[enabled_5_MR %in% common_tfs_enabled_list]
enabled_5_MR_and_TF <- paste0(enabled_5_MR_and_TF,collapse = ", ")
enabled_4_MR <- enabled_df[enabled_df$No_ICR_Cancers==4,]$MR
enabled_4_MR_and_TF <- enabled_4_MR[enabled_4_MR %in% common_tfs_enabled_list]
enabled_4_MR_and_TF <- paste0(enabled_4_MR_and_TF,collapse = ", ")
temp_enabled_df <- enabled_df[enabled_df$No_ICR_Cancers>=4,]
temp_enabled_df <- temp_enabled_df[temp_enabled_df$MR %in% common_tfs_enabled_list, ]
write.table(temp_enabled_df,"Results/Text_Results/All_MRs_TFs_ICR_Enabled_Supplementary_Table_S6.csv",row.names=F,col.names=T,quote=F,sep=",")


#Get list of common MRs which are also TFs in 4 out of 8 cancers
common_mrs_enabled_list <- temp_enabled_df$MR;
common_mrs_enabled_list <- sort(common_mrs_enabled_list)
common_mrs_enabled_activity_matrix <- matrix(0,nrow=length(common_mrs_enabled_list),ncol=2*length(icr_enabled))
rownames(common_mrs_enabled_activity_matrix) <- common_mrs_enabled_list
colnames(common_mrs_enabled_activity_matrix) <- c(icr_enabled,icr_enabled)
for (i in 1:length(icr_enabled))
{
  cancer <- icr_enabled[i]
  temp_df <- output_df[output_df$TF %in% common_mrs_enabled_list & output_df$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$TF),]
  icr_high <- temp_df$Avg_ICR_High
  icr_low <- temp_df$Avg_ICR_Low
  common_mrs_enabled_activity_matrix[common_mrs_enabled_list,i] <- icr_high
  common_mrs_enabled_activity_matrix[common_mrs_enabled_list,(i+length(icr_enabled))] <- icr_low
}

colcol <- matrix(0,nrow=2*length(icr_enabled),ncol=1)
colcol[c(1:length(icr_enabled)),1] <- "yellow"
colcol[c((length(icr_enabled)+1):(2*length(icr_enabled))),1] <- "green"

#Figure 4A 
#=========================================================================================================
pdf("Results/Figures/Common_TopMRs_Activity_ICR_Enabled_Figure_4A.pdf",height = 12, width=12, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.75)
p1 <- heatmap.3(common_mrs_enabled_activity_matrix, Rowv = TRUE, Colv=TRUE, col = bluered(100), scale="none", main= "Activity of Common MRs in ICR Enabled Cancers", # (>=4 out of 8)",
          dendrogram = "both", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.0, cexCol = 2.0, cellnote = ifelse(common_mrs_enabled_activity_matrix>0,"+","-"), notecex = 1, notecol = "black")
dev.off()
ordered_mrs_enabled <- rev(rownames(common_mrs_enabled_activity_matrix)[p1$rowInd])
nes_info <- NULL
for (i in 1:length(ordered_mrs_enabled))
{
  nes_info <- c(nes_info,mean(output_df[output_df$TF==ordered_mrs_enabled[i]  & output_df$Cancer %in% icr_enabled,]$NES))
}
names(nes_info) <- ordered_mrs_enabled

#Get activity of mrs of interest from enabled cancers (From "TGFB1I1" onwards)
mrs_of_interest_enabled <- names(nes_info[nes_info<0])
enabled_activity_df <- NULL
for (cancer in icr_enabled)
{
  load(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  
  high_indices_table <- read.table(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_low_indices.csv"),header=TRUE)
  high_indices <- high_indices_table$x
  low_indices <- low_indices_table$x  
  
  for (topmr in mrs_of_interest_enabled)
  {
    high_activity <- amat[topmr,high_indices]
    low_activity <- amat[topmr,low_indices]
    enabled_activity_df <- rbind(enabled_activity_df,
                                 cbind(rep(cancer,length(high_activity)+length(low_activity)),
                                       rep(topmr,length(high_activity)+length(low_activity)),
                                       c(rep("ICR High",length(high_activity)),rep("ICR Low",length(low_activity))),
                                       c(high_activity,low_activity)))
  }
}

enabled_activity_df <- as.data.frame(enabled_activity_df)
colnames(enabled_activity_df) <- c("Cancer","TopMR","Phenotype","Activity_Value")
enabled_activity_df$Cancer <- as.character(as.vector(enabled_activity_df$Cancer))
enabled_activity_df$TopMR <- as.character(as.vector(enabled_activity_df$TopMR))
enabled_activity_df$Phenotype <- as.character(as.vector(enabled_activity_df$Phenotype))
enabled_activity_df$Activity_Value <- as.numeric(as.vector(enabled_activity_df$Activity_Value))

supp_p3 <- ggplot(data = enabled_activity_df, aes(x=Cancer, y=Activity_Value)) + 
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ TopMR, nrow=5, ncol=5) + xlab("Cancer") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top MRs of interest in ICR Enabled Cancers (>=4 out of 8)") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))
#ggsave(file="Results/Figures/Supp_Figure3.jpg",plot = supp_p3, device = jpeg(), width = 14, height=12, units = "in", dpi = 300)
#dev.off()

mrs_of_interest_enabled_high <- names(nes_info[nes_info>0])
enabled_high_activity_df <- NULL
for (cancer in icr_enabled)
{
  load(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  
  high_indices_table <- read.table(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_low_indices.csv"),header=TRUE)
  high_indices <- high_indices_table$x
  low_indices <- low_indices_table$x  
  
  for (topmr in mrs_of_interest_enabled_high)
  {
    high_activity <- amat[topmr,high_indices]
    low_activity <- amat[topmr,low_indices]
    enabled_high_activity_df <- rbind(enabled_high_activity_df,
                                 cbind(rep(cancer,length(high_activity)+length(low_activity)),
                                       rep(topmr,length(high_activity)+length(low_activity)),
                                       c(rep("ICR High",length(high_activity)),rep("ICR Low",length(low_activity))),
                                       c(high_activity,low_activity)))
  }
}

enabled_high_activity_df <- as.data.frame(enabled_high_activity_df)
colnames(enabled_high_activity_df) <- c("Cancer","TopMR","Phenotype","Activity_Value")
enabled_high_activity_df$Cancer <- as.character(as.vector(enabled_high_activity_df$Cancer))
enabled_high_activity_df$TopMR <- as.character(as.vector(enabled_high_activity_df$TopMR))
enabled_high_activity_df$Phenotype <- as.character(as.vector(enabled_high_activity_df$Phenotype))
enabled_high_activity_df$Activity_Value <- as.numeric(as.vector(enabled_high_activity_df$Activity_Value))

supp_p3_b <- ggplot(data = enabled_high_activity_df, aes(x=Cancer, y=Activity_Value)) + 
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ TopMR, nrow=9, ncol=8) + xlab("Cancer") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top MRs of interest in ICR Enabled Cancers (>=4 out of 8)") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))
#ggsave(file="Results/Figures/Supp_Figure3.jpg",plot = supp_p3, device = jpeg(), width = 14, height=12, units = "in", dpi = 300)
#dev.off()

#==========================================================================================================================
###########################################################################################################################
common_tfs_disabled_list <- c(output_df[output_df$Cancer==icr_disabled[1],]$TF)
for (i in 2:length(icr_disabled))
{
  new_tfs <- output_df[output_df$Cancer==icr_disabled[i],]$TF
  common_tfs_disabled_list <- intersect(common_tfs_disabled_list,new_tfs)
}

first_cancer <- icr_disabled[1];
load(paste0("Results/",first_cancer,"/Adjacency_Matrix/",first_cancer,"_Full_TopMR_Info_FGSEA_BC.Rdata"))
topmr <- topmr_info$pathway
list_common_mrs_disabled <- as.character(as.vector(topmr));
topMRs_disabled <- rep(list("LGG"),length(list_common_mrs_disabled))
names(topMRs_disabled) <- list_common_mrs_disabled
for(i in 2:length(icr_disabled))
{
  cancer_type <- icr_disabled[i]
  load(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_TopMR_Info_FGSEA_BC.Rdata"))
  list_common_mrs_disabled <- as.character(as.vector(topmr_info$pathway));
  for (j in 1:length(list_common_mrs_disabled))
  {
    if (list_common_mrs_disabled[j] %in% names(topMRs_disabled))
    {
      topMRs_disabled[[list_common_mrs_disabled[j]]] <- c(topMRs_disabled[[list_common_mrs_disabled[j]]],cancer_type)
    } else {
      topMRs_disabled[[list_common_mrs_disabled[j]]] <- cancer_type
    }
  }
}

all_icr_cancers <- unlist(lapply(topMRs_disabled,function(x) paste(x,collapse=" ")))
no_icr_cancers <- unlist(lapply(topMRs_disabled,function(x) length(x)))
disabled_df <- cbind(names(all_icr_cancers),as.character(all_icr_cancers),as.numeric(no_icr_cancers))
disabled_df <- as.data.frame(disabled_df)
colnames(disabled_df) <- c("MR","List_ICR_Cancers","No_ICR_Cancers")
disabled_df$MR <- as.character(as.vector(disabled_df$MR))
disabled_df$List_ICR_Cancers <- as.character(as.vector(disabled_df$List_ICR_Cancers))
disabled_df$No_ICR_Cancers <- as.character(as.vector(disabled_df$No_ICR_Cancers))
disabled_df <- disabled_df[order(disabled_df$No_ICR_Cancers,decreasing=T),]

disabled_4_MR_and_TF <- disabled_df[disabled_df$No_ICR_Cancers==4,]$MR
disabled_3_MR <- disabled_df[disabled_df$No_ICR_Cancers==3,]$MR
disabled_3_MR_and_TF <- disabled_3_MR[disabled_3_MR %in% common_tfs_disabled_list]
disabled_3_MR_and_TF <- paste0(disabled_3_MR_and_TF,collapse = ", ")
disabled_2_MR <- disabled_df[disabled_df$No_ICR_Cancers==2,]$MR
disabled_2_MR_and_TF <- disabled_2_MR[disabled_2_MR %in% common_tfs_disabled_list]
disabled_2_MR_and_TF <- paste0(disabled_2_MR_and_TF,collapse = ", ")
disabled_1_MR <- disabled_df[disabled_df$No_ICR_Cancers==1,]$MR
disabled_1_MR_and_TF <- disabled_1_MR[disabled_1_MR %in% common_tfs_enabled_list]
disabled_1_MR_and_TF <- paste0(disabled_1_MR_and_TF,collapse = ", ")
temp_disabled_df <- disabled_df[disabled_df$No_ICR_Cancers>=2,]
temp_disabled_df <- temp_disabled_df[temp_disabled_df$MR %in% common_tfs_disabled_list, ]
write.table(temp_disabled_df,"Results/Text_Results/All_MRs_TFs_ICR_Enabled_Supplementary_Table_S7.csv",row.names=F,col.names=T,quote=F,sep=",")


#Get list of common MRs which are also TFs in 2 out of 4 cancers
common_mrs_disabled_list <- temp_disabled_df$MR;
common_mrs_disabled_list <- sort(common_mrs_disabled_list)
common_mrs_disabled_activity_matrix <- matrix(0,nrow=length(common_mrs_disabled_list),ncol=2*length(icr_disabled))
rownames(common_mrs_disabled_activity_matrix) <- common_mrs_disabled_list
colnames(common_mrs_disabled_activity_matrix) <- c(icr_disabled,icr_disabled)
for (i in 1:length(icr_disabled))
{
  cancer <- icr_disabled[i]
  temp_df <- output_df[output_df$TF %in% common_mrs_disabled_list & output_df$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$TF),]
  icr_high <- temp_df$Avg_ICR_High
  icr_low <- temp_df$Avg_ICR_Low
  common_mrs_disabled_activity_matrix[common_mrs_disabled_list,i] <- icr_high
  common_mrs_disabled_activity_matrix[common_mrs_disabled_list,(i+length(icr_disabled))] <- icr_low
}

colcol <- matrix(0,nrow=2*length(icr_disabled),ncol=1)
colcol[c(1:length(icr_disabled)),1] <- "yellow"
colcol[c((length(icr_disabled)+1):(2*length(icr_disabled))),1] <- "green"

#Figure 4B
#=========================================================================================================
pdf("Results/Figures/Common_TopMRs_Activity_ICR_Disabled_Figure_4B.pdf",height = 12, width=12, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.65)
p2 <- heatmap.3(common_mrs_disabled_activity_matrix, Rowv = TRUE, Colv=TRUE, col = bluered(100), scale="none", main= "Activity of Common MRs in ICR Disabled Cancers", # (>=2 out of 4)",
                dendrogram = "both", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.5, cexCol = 2, cellnote = ifelse(common_mrs_disabled_activity_matrix>0,"+","-"), notecex = 0.75, notecol = "black")
dev.off()

ordered_mrs_disabled <- rev(rownames(common_mrs_disabled_activity_matrix)[p2$rowInd])
nes_info_disabled <- NULL
for (i in 1:length(ordered_mrs_disabled))
{
  nes_info_disabled <- c(nes_info_disabled,mean(output_df[output_df$TF==ordered_mrs_disabled[i] & output_df$Cancer %in% icr_disabled,]$NES))
}
names(nes_info_disabled) <- ordered_mrs_disabled


#Get activity of mrs of interest from disabled cancers
mrs_of_interest_disabled <- names(nes_info_disabled[nes_info_disabled<0])
disabled_activity_df <- NULL
for (cancer in icr_disabled)
{
  load(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  
  high_indices_table <- read.table(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_low_indices.csv"),header=TRUE)
  high_indices <- high_indices_table$x
  low_indices <- low_indices_table$x  
  
  for (topmr in mrs_of_interest_disabled)
  {
    high_activity <- amat[topmr,high_indices]
    low_activity <- amat[topmr,low_indices]
    disabled_activity_df <- rbind(disabled_activity_df,
                                 cbind(rep(cancer,length(high_activity)+length(low_activity)),
                                       rep(topmr,length(high_activity)+length(low_activity)),
                                       c(rep("ICR High",length(high_activity)),rep("ICR Low",length(low_activity))),
                                       c(high_activity,low_activity)))
  }
}

disabled_activity_df <- as.data.frame(disabled_activity_df)
colnames(disabled_activity_df) <- c("Cancer","TopMR","Phenotype","Activity_Value")
disabled_activity_df$Cancer <- as.character(as.vector(disabled_activity_df$Cancer))
disabled_activity_df$TopMR <- as.character(as.vector(disabled_activity_df$TopMR))
disabled_activity_df$Phenotype <- as.character(as.vector(disabled_activity_df$Phenotype))
disabled_activity_df$Activity_Value <- as.numeric(as.vector(disabled_activity_df$Activity_Value))

supp_p4 <- ggplot(data = disabled_activity_df, aes(x=Cancer, y=Activity_Value)) + 
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ TopMR, nrow=9, ncol=9) + xlab("Cancer") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top MRs of interest in ICR Enabled Cancers (>=2 out of 4)") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))
#ggsave(file="Results/Figures/Supp_Figure4.jpg",plot = supp_p4, device = jpeg(), width = 14, height=12, units = "in", dpi = 300)
#dev.off()

#=========================================================================================================================
###########################################################################################################################
common_mrs_icr_high_all_cancers <- union(names(nes_info[nes_info>0]),names(nes_info_disabled[nes_info_disabled>0]))
common_mrs_icr_low_all_cancers <- union(names(nes_info[nes_info<0]),names(nes_info_disabled[nes_info_disabled<0]))
all_cancers <- c(icr_enabled,icr_disabled)
common_mrs_icr_high_activity_matrix <- matrix(0,nrow=length(common_mrs_icr_high_all_cancers),ncol=(length(icr_enabled)+length(icr_disabled)))
common_mrs_icr_high_nes_matrix <- matrix(0,nrow=length(common_mrs_icr_high_all_cancers),ncol=(length(icr_enabled)+length(icr_disabled)))
common_mrs_df <- NULL
for (i in 1:length(common_mrs_icr_high_all_cancers))
{
  topmr <- common_mrs_icr_high_all_cancers[i]
  for (j in 1:length(all_cancers))
  {
    cancer <- all_cancers[j]
    temp <- output_df[output_df$TF==topmr & output_df$Cancer==cancer,]
    if (dim(temp)[1]>0)
    {
      activity_val <- temp$Avg_ICR_High
      common_mrs_icr_high_activity_matrix[i,j] <- activity_val
      common_mrs_icr_high_nes_matrix[i,j] <- temp$NES
      temp_df <- cbind(cancer,topmr,temp$Padj,temp$NES,activity_val)
    }
    else{
      temp_df <- cbind(cancer,topmr,1,0,0)
    }
    common_mrs_df <- rbind(common_mrs_df,temp_df)
  }
  
}
rownames(common_mrs_icr_high_activity_matrix) <- common_mrs_icr_high_all_cancers
colnames(common_mrs_icr_high_activity_matrix) <- all_cancers
rownames(common_mrs_icr_high_nes_matrix) <- common_mrs_icr_high_all_cancers
colnames(common_mrs_icr_high_nes_matrix) <- all_cancers
common_mrs_df <- as.data.frame(common_mrs_df)
colnames(common_mrs_df) <- c("Cancer","TopMRs","Padj","NES","Avg_ICR_High_Activity")
common_mrs_df$Cancer <- as.character(as.vector(common_mrs_df$Cancer))
common_mrs_df$TopMRs <- as.character(as.vector(common_mrs_df$TopMRs))
common_mrs_df$Padj <- as.numeric(as.vector(common_mrs_df$Padj))
common_mrs_df$NES <- as.numeric(as.vector(common_mrs_df$NES))
common_mrs_df$Avg_ICR_High_Activity <- as.numeric(as.vector(common_mrs_df$Avg_ICR_High_Activity))

avg_icr_high_nes_all_cancers <- rowMeans(common_mrs_icr_high_nes_matrix)
indices <- order(avg_icr_high_nes_all_cancers,decreasing=TRUE)
common_mrs_icr_high_activity_matrix <- common_mrs_icr_high_activity_matrix[indices,]
common_mrs_icr_high_nes_matrix <- common_mrs_icr_high_nes_matrix[indices,]

#Supplementary Figure S5
#==============================================================
icr_high_df <- cbind(common_mrs_icr_high_activity_matrix,common_mrs_icr_high_nes_matrix)
icr_high_df <- as.data.frame(icr_high_df)
colnames(icr_high_df) <- c(paste0(all_cancers,"_Avg_ICR_High_Activity"),paste0(all_cancers,"_Avg_NES"))
avg_nes_vec <- rowMeans(icr_high_df[,c((length(all_cancers)+1):ncol(icr_high_df))])
nes_more_than_0 <- rowSums(icr_high_df[,c((length(all_cancers)+1):ncol(icr_high_df))]>0)
nes_less_than_0 <- rowSums(icr_high_df[,c((length(all_cancers)+1):ncol(icr_high_df))]<0)
nes_equal_to_0 <- rowSums(icr_high_df[,c((length(all_cancers)+1):ncol(icr_high_df))]==0)
icr_high_df$"Avg_NES_All_Cancers" <- avg_nes_vec
icr_high_df$"NES>0" <- nes_more_than_0
icr_high_df$"NES<0" <- nes_less_than_0
icr_high_df$"NES=0" <- nes_equal_to_0

diff_mrs_icr_high <- names(which(rowSums(icr_high_df[,c(1:8)]<0)>4 & rowSums(icr_high_df[,c(9:12)]>0)>2))
common_mrs_icr_high_final <- setdiff(rownames(common_mrs_icr_high_nes_matrix[rowMeans(common_mrs_icr_high_nes_matrix)>0.1,]),diff_mrs_icr_high)
write.table(icr_high_df[common_mrs_icr_high_final,],"Results/Text_Results/All_MRS_ICR_High_Supplementary_Table_S6.csv",row.names=T,col.names=T,quote=F,sep=",")


colcol <- matrix(0,nrow=length(all_cancers),ncol=1)
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_cancers)),1] <- "#660066"

#Figure 4C
#====================================================================================
pdf("Results/Figures/Common_TopMRs_Activity_ICR_High_Figure_4C.pdf",height = 10, width=12, pointsize = 12)
p3 <- heatmap.3(common_mrs_icr_high_activity_matrix[common_mrs_icr_high_final,], Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Activity of MRs specific to ICR High Phenotype for 12 cancers in ICR High samples",
                dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.5, cexCol = 1.5, cellnote = ifelse(common_mrs_icr_high_activity_matrix[common_mrs_icr_high_final,]>=0,"+","-"), notecex = 1.0, notecol = "black")
dev.off()


diff_mrs_icr_high_activity_matrix <- common_mrs_icr_high_activity_matrix[diff_mrs_icr_high,]
diff_mrs_icr_high_pval_matrix <- matrix(1,nrow=nrow(diff_mrs_icr_high_activity_matrix),ncol=ncol(diff_mrs_icr_high_activity_matrix))
rownames(diff_mrs_icr_high_pval_matrix) <- rownames(diff_mrs_icr_high_activity_matrix)
colnames(diff_mrs_icr_high_pval_matrix) <- colnames(diff_mrs_icr_high_activity_matrix)
for (i in 1:nrow(diff_mrs_icr_high_pval_matrix))
{
  mr <- rownames(diff_mrs_icr_high_pval_matrix)[i]
  for (j in 1:ncol(diff_mrs_icr_high_pval_matrix))
  {
    cancer <- colnames(diff_mrs_icr_high_pval_matrix)[j]
    if (nrow(output_df[output_df$TF==mr & output_df$Cancer==cancer,])>0)
    {
      diff_mrs_icr_high_pval_matrix[i,j] <- output_df[output_df$TF==mr & output_df$Cancer==cancer,]$Padj
    }
  }
}
hc.rows <- hclust(dist(diff_mrs_icr_high_activity_matrix),method='ward.D2')
colcol <- matrix(0,nrow=length(all_cancers),ncol=2)
colnames(colcol) <- c("Enabled/Disabled Cancers","ICR High")
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_cancers)),1] <- "#660066"
colcol[c(1:12),2] <- "yellow"

pdf("Results/Figures/Different_TopMRs_Activity_ICR_High.pdf",height = 10, width=12, pointsize = 12)
par(bg="black")
par(fg="white",col.axis="white",col.main="white",col.lab="white", cex.main=1.25)
heatmap.3(diff_mrs_icr_high_activity_matrix[hc.rows$order,], Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Activity of ICR High MRs different between ICR Enabled & ICR Disabled cancers",
                dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 3,
                margins = c(6,6), useRaster = FALSE, cexRow = 1.5, cexCol = 1.5, cellnote = ifelse(diff_mrs_icr_high_pval_matrix[hc.rows$order,]<=0.05,"*",""), notecex = 2, notecol = "black")
dev.off()

#==================================================================================================================================
#==============================================================================================================================
common_mrs_icr_low_all_cancers <- union(names(nes_info[nes_info<0]),names(nes_info_disabled[nes_info_disabled<0]))
all_cancers <- c(icr_enabled,icr_disabled)
common_mrs_icr_low_activity_matrix <- matrix(0,nrow=length(common_mrs_icr_low_all_cancers),ncol=(length(icr_enabled)+length(icr_disabled)))
common_mrs_icr_low_nes_matrix <- matrix(0,nrow=length(common_mrs_icr_low_all_cancers),ncol=(length(icr_enabled)+length(icr_disabled)))
for (i in 1:length(common_mrs_icr_low_all_cancers))
{
  topmr <- common_mrs_icr_low_all_cancers[i]
  for (j in 1:length(all_cancers))
  {
    cancer <- all_cancers[j]
    temp <- output_df[output_df$TF==topmr & output_df$Cancer==cancer,]
    if (dim(temp)[1]>0)
    {
      activity_val <- temp$Avg_ICR_Low
      common_mrs_icr_low_activity_matrix[i,j] <- activity_val
      common_mrs_icr_low_nes_matrix[i,j] <- temp$NES
    }
  }
}
rownames(common_mrs_icr_low_activity_matrix) <- common_mrs_icr_low_all_cancers
colnames(common_mrs_icr_low_activity_matrix) <- all_cancers
rownames(common_mrs_icr_low_nes_matrix) <- common_mrs_icr_low_all_cancers
colnames(common_mrs_icr_low_nes_matrix) <- all_cancers
avg_icr_low_nes_all_cancers <- rowMeans(common_mrs_icr_low_nes_matrix)
indices <- order(avg_icr_low_nes_all_cancers,decreasing=FALSE)
common_mrs_icr_low_activity_matrix <- common_mrs_icr_low_activity_matrix[indices,]
common_mrs_icr_low_nes_matrix <- common_mrs_icr_low_nes_matrix[indices,]

#Supplementary Figure S6
#===============================
icr_low_df <- cbind(common_mrs_icr_low_activity_matrix,common_mrs_icr_low_nes_matrix)
icr_low_df <- as.data.frame(icr_low_df)
colnames(icr_low_df) <- c(paste0(all_cancers,"_Avg_ICR_Low_Activity"),paste0(all_cancers,"_Avg_NES"))
avg_nes_vec <- rowMeans(icr_low_df[,c((length(all_cancers)+1):ncol(icr_low_df))])
nes_more_than_0 <- rowSums(icr_low_df[,c((length(all_cancers)+1):ncol(icr_low_df))]>0)
nes_less_than_0 <- rowSums(icr_low_df[,c((length(all_cancers)+1):ncol(icr_low_df))]<0)
nes_equal_to_0 <- rowSums(icr_low_df[,c((length(all_cancers)+1):ncol(icr_low_df))]==0)
icr_low_df$"Avg_NES_All_Cancers" <- avg_nes_vec
icr_low_df$"NES>0" <- nes_more_than_0
icr_low_df$"NES<0" <- nes_less_than_0
icr_low_df$"NES=0" <- nes_equal_to_0

diff_mrs_icr_low <- names(which(rowSums(icr_low_df[,c(1:8)]>0)>4 & rowSums(icr_low_df[,c(9:12)]<0)>2))
common_mrs_icr_low_final <- setdiff(rownames(common_mrs_icr_low_nes_matrix[rowMeans(common_mrs_icr_low_nes_matrix)< -0.1,]),diff_mrs_icr_low)
write.table(icr_low_df[common_mrs_icr_low_final,],"Results/Text_Results/All_MRS_ICR_Low_Supplementary_Table_S7.csv",row.names=T,col.names=T,quote=F,sep=",")

colcol <- matrix(0,nrow=length(all_cancers),ncol=1)
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_cancers)),1] <- "#660066"

#Figure 4D
#====================================================================================
pdf("Results/Figures/Common_TopMRs_Activity_ICR_Low_Figure_4D.pdf",height = 10, width=12, pointsize = 12)
p4 <- heatmap.3(common_mrs_icr_low_activity_matrix[common_mrs_icr_low_final,], Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Activity of MRs specific to ICR Low Phenotype for 12 cancers in ICR Low samples",
                dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.75, cexCol = 1.5, cellnote = ifelse(common_mrs_icr_low_activity_matrix[common_mrs_icr_low_final,]>=0,"+","-"), notecex = 1.0, notecol = "black")
dev.off()


diff_mrs_icr_low_activity_matrix <- common_mrs_icr_low_activity_matrix[diff_mrs_icr_low,]
diff_mrs_icr_low_pval_matrix <- matrix(1,nrow=nrow(diff_mrs_icr_low_activity_matrix),ncol=ncol(diff_mrs_icr_low_activity_matrix))
rownames(diff_mrs_icr_low_pval_matrix) <- rownames(diff_mrs_icr_low_activity_matrix)
colnames(diff_mrs_icr_low_pval_matrix) <- colnames(diff_mrs_icr_low_activity_matrix)
for (i in 1:nrow(diff_mrs_icr_low_pval_matrix))
{
  mr <- rownames(diff_mrs_icr_low_pval_matrix)[i]
  for (j in 1:ncol(diff_mrs_icr_low_pval_matrix))
  {
    cancer <- colnames(diff_mrs_icr_low_pval_matrix)[j]
    if (nrow(output_df[output_df$TF==mr & output_df$Cancer==cancer,])>0)
    {
      diff_mrs_icr_low_pval_matrix[i,j] <- output_df[output_df$TF==mr & output_df$Cancer==cancer,]$Padj
    }
  }
}
hc.rows <- hclust(dist(diff_mrs_icr_low_activity_matrix),method='ward.D2')
colcol <- matrix(0,nrow=length(all_cancers),ncol=2)
colnames(colcol) <- c("Enabled/Disabled Cancers","ICR Low")
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_cancers)),1] <- "#660066"
colcol[c(1:12),2] <- "green"

pdf("Results/Figures/Different_TopMRs_Activity_ICR_Low.pdf",height = 10, width=12, pointsize = 12)
par(bg="black")
par(fg="white",col.axis="white",col.main="white",col.lab="white", cex.main=1.25)
heatmap.3(diff_mrs_icr_low_activity_matrix[hc.rows$order,], Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Activity of ICR Low MRs different between ICR Enabled & ICR Disabled cancers",
          dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 3,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.5, cexCol = 1.5, cellnote = ifelse(diff_mrs_icr_low_pval_matrix[hc.rows$order,]<=0.05,"*",""), notecex = 2, notecol = "black")
dev.off()

##################################################################################################################################

# #Gene list of interest for ICR Low phenotype
# entrez_icr_low <- mapIds(org.Hs.eg.db, common_mrs_icr_low_final, 'ENTREZID', 'SYMBOL')
# entrez_all <- mapIds(org.Hs.eg.db,target_gene_names,'ENTREZID','SYMBOL')
# wp2gene <- read.gmt("Results/Text_Results/wikipathways-20190810-gmt-Homo_sapiens.gmt")
# wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
# wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
# wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
# m_df <- msigdbr(species = "Homo sapiens")
# m_t2g <- msigdbr(species = "Homo sapiens", category="C7") %>% 
#   dplyr::select(gs_name, entrez_gene)
# 
# library(enrichplot)
# ewp <- enricher(entrez_icr_low, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
# em <- enricher(entrez_icr_low, TERM2GENE=m_t2g)

icr_high_goterms_df <- read.table("Results/Text_Results/Enriched_GOTerms_ICR_High.csv",header=TRUE,sep="\t")
temp_df <- icr_high_goterms_df[-log10(icr_high_goterms_df$p.value)>40,]

supp_p7a <- ggplot(icr_high_goterms_df,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) + geom_point(aes(color=as.factor(term_level))) + 
           xlab("Percentage of MRs in each GO Term") + ylab("-log10(Pvalue)") + scale_shape_manual(name = "GO Terms", values = 1:nlevels(as.factor(icr_high_goterms_df$term_category))) +
            scale_color_manual(name = "Term Level", values =  c("grey","violet","orange")) +
            geom_hline(yintercept=c(0,20,40), color=c("black","blue","red")) + geom_vline(xintercept = 0) + 
            theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlim(c(0.0,70))+
            geom_text(aes(x=100*generatio,y=-log10(p.value),label=term_name),
              data=temp_df,
              hjust = 0, vjust = 0, nudge_x = 0.025, angle = 0,
              size = 3, check_overlap = T) +
            ggtitle("Enriched GO Terms for ICR High Phenotype") + theme(text = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="Results/Figures/Supp_Figure7A.jpg",plot = supp_p7a, device = jpeg(), width = 7, height=6, units = "in", dpi = 300)
dev.off()

icr_low_goterms_df <- read.table("Results/Text_Results/Enriched_GOTerms_ICR_Low.csv",header=TRUE,sep="\t")
temp_df <- icr_low_goterms_df[-log10(icr_low_goterms_df$p.value)>38,]

supp_p7b <- ggplot(icr_low_goterms_df,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) + geom_point(aes(color=as.factor(term_level))) + 
  xlab("Percentage of MRs in each GO Term") + ylab("-log10(Pvalue)") + scale_shape_manual(name = "GO Terms", values = 1:nlevels(as.factor(icr_low_goterms_df$term_category))) +
  scale_color_manual(name = "Term Level", values =  c("grey","violet","orange")) +
  geom_hline(yintercept=c(0,20,40), color=c("black","blue","red")) + geom_vline(xintercept = 0) + 
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(aes(x=100*generatio,y=-log10(p.value),label=term_name),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.025, angle = 0,
            size = 3, check_overlap = T) +
  coord_cartesian(xlim =c(0, 15), ylim = c(0, 50)) +
  ggtitle("Enriched GO Terms for ICR Low Phenotype") + theme(text = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="Results/Figures/Supp_Figure7B.jpg",plot = supp_p7b, device = jpeg(), width = 7, height=6, units = "in", dpi = 300)
dev.off()

#Make Figure 5A_2
#===========================================================================================================
colors <- c('#3cb44b', '#4363d8', '#f58231', '#911eb4', '#f032e6', '#46f0f0' , '#bcf60c', '#a9a9a9', '#e6194b', '#ffe119', '#ffe119', '#e6beff')
icr_low_pathways_df <- read.table("Results/Text_Results/Enriched_Pathways_ICR_Low.csv",header=TRUE,sep="\t")
icr_low_pathways_df$Description <- as.character(as.vector(icr_low_pathways_df$Description))
icr_low_pathways_df$GeneRatio <- as.character(as.vector(icr_low_pathways_df$GeneRatio))
icr_low_pathways_df$genes <- as.character(as.vector(icr_low_pathways_df$genes))
icr_low_pathways_df$Description <- paste0(icr_low_pathways_df$Description," [",icr_low_pathways_df$GeneRatio,"] ")
icr_low_pathways_df <- icr_low_pathways_df[order(icr_low_pathways_df$generatio),]
p5 <- ggplot(data=icr_low_pathways_df,aes(x=generatio,y=reorder(Description,generatio),size=-log10(pvalue)))+
      geom_point(aes(color=-log10(pvalue)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
      scale_size_continuous(name = "-log10(P.adjust)")+   guides(color=guide_legend(), size = guide_legend())+
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background=element_rect(fill = "black"),axis.line = element_line(colour = "white"),
      #  panel.background = element_rect(fill = 'black'),legend.background = element_rect(fill = "black", color = NA),
      #  legend.key = element_rect(color = "gray", fill = "black"), legend.title = element_text(color = "white"), legend.text = element_text(color = "white"))+
      ggtitle("Enriched Pathways for ICR Low Phenotype") + theme(text = element_text(size=16,color="black")) + theme(axis.text = element_text(size=14, color="black")) + theme(plot.title = element_text(hjust = 0.5,color="black")) +
      theme(axis.text.y=element_text(color=colors[(icr_low_pathways_df$clusters)]))
ggsave(filename = "Results/Figures/Enriched_Pathways_ICR_Low_Figure5A_2.pdf",plot=p5,device=pdf(), height=12, width=15, units="in", dpi=300)
dev.off()
  
#new_col=c("red","green","blue")
#Supplementary Figure S8.A
#================================================================================================================
icr_high_pathways_df <- read.table("Results/Text_Results/Enriched_Pathways_ICR_High.csv",header=TRUE,sep="\t")
icr_high_pathways_df$Description <- as.character(as.vector(icr_high_pathways_df$Description))
icr_high_pathways_df$GeneRatio <- as.character(as.vector(icr_high_pathways_df$GeneRatio))
icr_high_pathways_df$Description <- paste0(icr_high_pathways_df$Description," [",icr_high_pathways_df$GeneRatio,"] ")
icr_high_pathways_df$genes <- as.character(as.vector(icr_high_pathways_df$genes))
icr_high_pathways_df <- icr_high_pathways_df[order(icr_high_pathways_df$generatio),]
supp_p8a <- ggplot(data=icr_high_pathways_df,aes(x=generatio,y=reorder(Description,generatio),size=-log10(pvalue)))+
  geom_point(aes(color=-log10(pvalue)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(name = "-log10(P.adjust)")+ guides(color=guide_legend(), size = guide_legend())+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Enriched Pathways for ICR High Phenotype") + theme(text = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y=element_text(color=colors[as.factor(icr_high_pathways_df$clusters)]))
ggsave(filename = "Results/Figures/Supp_Figure8A.jpg",plot=supp_p8a,device=jpeg(), height=15, width=12, units="in", dpi=300)
dev.off()

#Supplementary Figure S8B
#=====================================================================================================================
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))

icr_high_pathways_df <- icr_high_pathways_df[order(icr_high_pathways_df$generatio,decreasing = T),]
icr_high_pathways_gene_involved <- matrix(0,nrow=nrow(icr_high_pathways_df),ncol=length(rownames(icr_high_df[which(icr_high_df$Avg_NES_All_Cancers>0),])))
rownames(icr_high_pathways_gene_involved) <- icr_high_pathways_df$Description
colnames(icr_high_pathways_gene_involved) <- rownames(icr_high_df[which(icr_high_df$Avg_NES_All_Cancers>0),])

for (i in 1:nrow(icr_high_pathways_df))
{
  genes_involved <- unlist(strsplit(icr_high_pathways_df[i,]$genes,split="; "))
  nes_scores <- icr_high_df[genes_involved,]$Avg_NES_All_Cancers
  icr_high_pathways_gene_involved[i,genes_involved] <- nes_scores
}

jpeg("Results/Figures/Supp_Figure8B.jpg",height = 900, width=1500, units="px",pointsize = 12)
Heatmap(icr_high_pathways_gene_involved,cluster_rows = FALSE, cluster_columns = FALSE, col=col_fun, row_names_side = "left", 
        name="NES", row_names_gp = gpar(fontsize = 6, col=colors[as.factor(icr_high_pathways_df$clusters)]), column_names_gp = gpar(fontsize = 8), column_title = "MRs in Enriched Pathways specific to ICR High")
dev.off()


#Figure 5A_1
#================================================================================================================
icr_low_pathways_df <- icr_low_pathways_df[order(icr_low_pathways_df$generatio,decreasing = T),]
icr_low_pathways_gene_involved <- matrix(0,nrow=nrow(icr_low_pathways_df),ncol=length(common_mrs_icr_low_final))
rownames(icr_low_pathways_gene_involved) <- icr_low_pathways_df$Description
colnames(icr_low_pathways_gene_involved) <- common_mrs_icr_low_final
for (i in 1:nrow(icr_low_pathways_df))
{
  genes_involved <- unlist(strsplit(icr_low_pathways_df[i,]$genes,split="; "))
  nes_scores <- icr_low_df[genes_involved,]$Avg_NES_All_Cancers
  icr_low_pathways_gene_involved[i,genes_involved] <- nes_scores
}

#Make the Sankey plot
icr_low_pathways_gene_involved_to_consider <- icr_low_pathways_gene_involved[,colSums(icr_low_pathways_gene_involved)<0]
ICR_Low <- list()
ICR_Low$nodes <- data.frame(name = c(colnames(icr_low_pathways_gene_involved_to_consider),rownames(icr_low_pathways_gene_involved_to_consider)))
nodesgroup <- c(rep("white",length(colnames(icr_low_pathways_gene_involved_to_consider))),colors[icr_low_pathways_df$clusters])
ICR_Low$nodes$nodesgroup <- nodesgroup
edgelist <- NULL
for (i in 1:nrow(icr_low_pathways_gene_involved_to_consider))
{
  for (j in 1:ncol(icr_low_pathways_gene_involved_to_consider))
  {
    pathway <- rownames(icr_low_pathways_gene_involved_to_consider)[i]
    mr <- colnames(icr_low_pathways_gene_involved_to_consider)[j]
    if (icr_low_pathways_gene_involved_to_consider[pathway,mr]<0)
    {
      mr_id <- which(ICR_Low$nodes$name==mr)-1
      pathway_id <- which(ICR_Low$nodes$name==pathway)-1
      temp <- cbind(mr_id,pathway_id,abs(icr_low_pathways_gene_involved_to_consider[pathway,mr]),ICR_Low$nodes[ICR_Low$nodes$name==pathway,"nodesgroup"])
      edgelist <- rbind(edgelist,temp)
    }
  }
}
edgelist <- as.data.frame(edgelist)
colnames(edgelist) <- c("source","target","value","type")
edgelist$source <- as.numeric(as.vector(edgelist$source))
edgelist$target <- as.numeric(as.vector(edgelist$target))
edgelist$value <- as.numeric(as.vector(edgelist$value))
edgelist$type <- as.factor(as.vector(edgelist$type))
ICR_Low$links <- edgelist

# putting in a data.frame might help see problems
color_scale <- data.frame(
  range = c(rep("white",length(colnames(icr_low_pathways_gene_involved_to_consider))),colors[icr_low_pathways_df$clusters]),
  domain = ICR_Low$nodes$name,
  nodes = ICR_Low$nodes,
  stringsAsFactors = FALSE
)

p <- sankeyNetwork(Links = ICR_Low$links, Nodes = ICR_Low$nodes, Source = "source",
                   Target = "target", Value = "value", LinkGroup = "type", NodeID = "name", NodeGroup = "nodesgroup",
                   units = "", fontSize = 16, nodeWidth = 25, iterations=0, fontFamily = "Arial" ) 


# p <- plot_ly(
#   type = "sankey",
#   domain = list(
#     x =  c(0,1),
#     y =  c(0,1)
#   ),
#   orientation = "h",
#   valueformat = ".0f",
#   valuesuffix = "",
#   
#   node = list(
#     label = ICR_Low$nodes$name,
#     color = ICR_Low$nodes$nodesgroup,
#     pad = 12,
#     thickness = 12,
#     line = list(
#       color = "black",
#       width = 0.5
#     )
#   ),
#   
#   link = list(
#     source = ICR_Low$links$source,
#     target = ICR_Low$links$target,
#     value =  ICR_Low$links$value,
#     color =  ICR_Low$links$type
#   )
#   
# )  %>% 
#   layout(
#     title = "Enriched Pathways for ICR Low Phenotype",
#     font = list(
#       size = 14,
#       color = "white"
#     ),
#     xaxis = list(showgrid = F, zeroline = F),
#     yaxis = list(showgrid = F, zeroline = F),
#     plot_bgcolor = 'blue',
#     paper_bgcolor = 'blue'
#   )
# 
# plotly_IMAGE(p, format = "png", file = "Results/Figures/Pdfs/Sankey_ICR_Low.png", width = 600, height = 480)




#Load the images 4A and 4B
im1 <- readPicture("Results/Figures/Common_TopMRs_Activity_ICR_Enabled_Figure_4A.svg")
p1 <- pictureGrob(im1)

im2 <- readPicture("Results/Figures/Common_TopMRs_Activity_ICR_Disabled_Figure_4B.svg")
p2 <- pictureGrob(im2)

im3 <- readPicture("Results/Figures/Common_TopMRs_Activity_ICR_High_Figure_4C.svg")
p3 <- pictureGrob(im3)

im4 <- readPicture("Results/Figures/Common_TopMRs_Activity_ICR_Low_Figure_4D.svg")
p4 <- pictureGrob(im4)

im5 <- readPicture("Results/Figures/Different_TopMRs_Activity_ICR_High_Figure_4E.svg")
p5 <- pictureGrob(im5)

im6 <- readPicture("Results/Figures/Different_TopMRs_Activity_ICR_Low_Figure_4F.svg")
p6 <- pictureGrob(im6)

g_final <- ggarrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,labels=c("A","B","C","D","E","F"),widths=c(2,2))
ggsave(filename = paste0("Results/Figures/Figure4.pdf"),device = pdf(), plot = g_final, dpi = 300, width = 16, height= 14, units = "in" )
dev.off()



