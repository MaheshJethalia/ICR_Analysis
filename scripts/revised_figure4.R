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
library(ComplexHeatmap)
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
library(networkD3)
library(matrixStats)
library(rbokeh)
warnings("off")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

registerDoMC(20)

#setwd('/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Old_Files/ICR_All_Info/scripts/')
setwd('.')

source('gene-reverse-network.R')
source('get_functions.R')


get_tf_median_activity <- function(cancer)
{
  load(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  high_indices_table <- read.table(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_low_indices.csv"),header=TRUE)
  high_indices <- high_indices_table$x
  low_indices <- low_indices_table$x 
  tfs_median_high_activity <- rowMedians(amat[,high_indices])
  tfs_median_low_activity <- rowMedians(amat[,low_indices])
  activity_info <- cbind(rownames(amat),rep(cancer,nrow(amat)),tfs_median_high_activity,tfs_median_low_activity)
  return(activity_info)
}

get_tf_activity_high_low <- function(cancer,tf)
{
  load(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  high_indices_table <- read.table(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_low_indices.csv"),header=TRUE)
  high_indices <- high_indices_table$x
  low_indices <- low_indices_table$x
  if (tf %in% rownames(amat))
  {
    tf_high_activity <- as.numeric(as.vector(amat[tf,high_indices]))
    tf_low_activity <- as.numeric(as.vector(amat[tf,low_indices]))
  }
  else{
    tf_high_activity <- rep(0,length(high_indices))
    tf_low_activity <- rep(0,length(low_indices))
  }
  return(list(tf_high_activity,tf_low_activity))
}


get_common_mrs <- function(cancer)
{
  load(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_TopMR_Info_FGSEA_BC_NES_1.Rdata"))
  load(paste0("../Results/",cancer,"/Combined/",cancer,"_Common_TopMRs.Rdata"))
  common_topmr <- topmr_info[topmr_info$pathway %in% common_top_mrs,]$pathway
  list_common_mrs <- as.character(as.vector(common_topmr));
  return(list_common_mrs)
}

icr_enabled <- c("BLCA","BRCA","HNSC","LIHC","SARC","SKCM","STAD","UCEC")
icr_disabled <- c("LGG","KIRC","PAAD","UVM")

tfs_enabled_activity_info <- NULL
tfs_disabled_activity_info <- NULL

#Get the list of TFs which are common to all the 12 ICR Enabled cancers and their median activity in ICR High vs ICR Low samples per cancer
#############################################################################################################################
colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#a9a9a9', '#008080', '#e6beff')
first_cancer <- icr_enabled[1];
tfs_enabled_activity_info <- get_tf_median_activity(first_cancer)
for (i in 2:length(icr_enabled))
{
  cancer <- icr_enabled[i]
  temp <- get_tf_median_activity(cancer)
  tfs_enabled_activity_info <- rbind(tfs_enabled_activity_info,temp)
}
tfs_enabled_activity_info <- as.data.frame(tfs_enabled_activity_info)
colnames(tfs_enabled_activity_info) <- c("TFs","Cancer","Median_ICR_High","Median_ICR_Low")
tfs_enabled_activity_info$TFs <- as.character(as.vector(tfs_enabled_activity_info$TFs))
tfs_enabled_activity_info$Cancer <- as.character(as.vector(tfs_enabled_activity_info$Cancer))
tfs_enabled_activity_info$Median_ICR_High <- as.numeric(as.vector(tfs_enabled_activity_info$Median_ICR_High))
tfs_enabled_activity_info$Median_ICR_Low <- as.numeric(as.vector(tfs_enabled_activity_info$Median_ICR_Low))

#Select those transcription regulators which are TFs with regulons of size >= 10 in at least 50% of the ICR Enabled Cancers
common_tfs_enabled_list <- names(which(table(tfs_enabled_activity_info$TFs)>7))
common_tfs_enabled_activity_info <- tfs_enabled_activity_info[tfs_enabled_activity_info$TFs %in% common_tfs_enabled_list,]

#Make a list of how many times does a common MR appear for the ICR Enabled cancers
first_cancer <- icr_enabled[1];
list_common_mrs_enabled <- get_common_mrs(first_cancer)
topMRs_enabled <- rep(list("BLCA"),length(list_common_mrs_enabled))
names(topMRs_enabled) <- list_common_mrs_enabled
for(i in 2:length(icr_enabled))
{
  cancer_type <- icr_enabled[i]
  list_common_mrs_enabled <- get_common_mrs(cancer_type)
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

enabled_7_MR <- enabled_df[enabled_df$No_ICR_Cancers==7,]$MR
#enabled_7_MR_and_TF <- enabled_7_MR[enabled_7_MR %in% common_tfs_enabled_list]
#enabled_7_MR_and_TF <- paste0(enabled_7_MR_and_TF,collapse = ", ")
temp_enabled_df <- enabled_df[enabled_df$No_ICR_Cancers>=4,]
temp_enabled_df <- temp_enabled_df[temp_enabled_df$MR %in% common_tfs_enabled_list, ]
write.table(temp_enabled_df,"../Results/Revised_Text_Results/All_MRs_TFs_ICR_Enabled_Supplementary_Table_Revised(v2)_S6.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get list of common MRs which are all TFs in 4 out of 8 cancers
#########################################################################################################################
common_mrs_enabled_list <- temp_enabled_df$MR;
common_mrs_enabled_list <- sort(common_mrs_enabled_list)
common_mrs_enabled_activity_matrix <- matrix(0,nrow=length(common_mrs_enabled_list),ncol=2*length(icr_enabled))
rownames(common_mrs_enabled_activity_matrix) <- common_mrs_enabled_list
colnames(common_mrs_enabled_activity_matrix) <- c(icr_enabled,icr_enabled)
for (i in 1:length(icr_enabled))
{
  cancer <- icr_enabled[i]
  temp_df <- common_tfs_enabled_activity_info[common_tfs_enabled_activity_info$TFs %in% common_mrs_enabled_list 
                                              & common_tfs_enabled_activity_info$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$TFs),]
  icr_high <- temp_df$Median_ICR_High
  icr_low <- temp_df$Median_ICR_Low
  common_mrs_enabled_activity_matrix[common_mrs_enabled_list,i] <- icr_high
  common_mrs_enabled_activity_matrix[common_mrs_enabled_list,(i+length(icr_enabled))] <- icr_low
}

colcol <- matrix(0,nrow=2*length(icr_enabled),ncol=1)
colcol[c(1:length(icr_enabled)),1] <- "yellow"
colcol[c((length(icr_enabled)+1):(2*length(icr_enabled))),1] <- "green"

#Perform Wilcox ranksum test or Mann-Whitney test to identify the MRs whose activity between the ICR High and ICR Low are significant for ICR Enabled cancers
#######################################################################################################
output_df <- common_tfs_enabled_activity_info
nes_info <- NULL
nes_pval <- NULL
for (i in 1:length(common_mrs_enabled_list))
{
  mr <- common_mrs_enabled_list[i]
  mr_icr_high_values <- NULL
  mr_icr_low_values <- NULL
  for (cancer in icr_enabled)
  {
      mr_activity_list <- get_tf_activity_high_low(cancer,mr)
      mr_icr_high_values <- c(mr_icr_high_values,mr_activity_list[[1]])
      mr_icr_low_values <- c(mr_icr_low_values,mr_activity_list[[2]])
  }
  nes_info <- c(nes_info,median(mr_icr_high_values)-median(mr_icr_low_values))
  nes_pval <- c(nes_pval,wilcox.test(mr_icr_high_values,mr_icr_low_values,exact=F)$p.value)
}
icr_enabled_median_comparison <- cbind(common_mrs_enabled_list,nes_info,nes_pval)
icr_enabled_median_comparison <- as.data.frame(icr_enabled_median_comparison)
colnames(icr_enabled_median_comparison) <- c("MR","FC_Median","Pval")
icr_enabled_median_comparison$MR <- as.character(as.vector(icr_enabled_median_comparison$MR))
icr_enabled_median_comparison$FC_Median <- round(as.numeric(as.vector(icr_enabled_median_comparison$FC_Median)),3)
icr_enabled_median_comparison$Pval <- p.adjust(as.numeric(as.vector(icr_enabled_median_comparison$Pval)),method = "fdr")
final_common_mrs_enabled_list <- icr_enabled_median_comparison[icr_enabled_median_comparison$Pval<0.05,]$MR

final_common_mrs_enabled_activity_matrix <- common_mrs_enabled_activity_matrix[final_common_mrs_enabled_list,]
#Figure 4A 
#=========================================================================================================
pdf("../Results/Revised_Figures/Pdfs/Common_TopMRs_Activity_ICR_Enabled_Figure_Revised(v2)_4A.pdf",height = 12, width=14, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.75)
p1 <- heatmap.3(final_common_mrs_enabled_activity_matrix, Rowv = TRUE, Colv=, col = bluered(100), scale="none", main= "Median Activity of Common MRs in ICR Enabled Cancers", # (>=4 out of 8)",
                dendrogram = "both", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.6, cexCol = 2.0, cellnote = ifelse(final_common_mrs_enabled_activity_matrix>0,"+","-"), notecex = 1, notecol = "black")
dev.off()

#Identify the MRs which are specific to ICR Low
ordered_mrs_enabled <- rev(rownames(final_common_mrs_enabled_activity_matrix)[p1$rowInd])
ordered_p_values_icr_enabled <- NULL
for (i in 1:length(ordered_mrs_enabled))
{
  mr <- ordered_mrs_enabled[i]
  adj_pval <- icr_enabled_median_comparison[icr_enabled_median_comparison$MR==mr,]$Pval
  temp <- cbind(mr,adj_pval)
  ordered_p_values_icr_enabled <- rbind(ordered_p_values_icr_enabled,temp)
}
ordered_p_values_icr_enabled <- as.data.frame(ordered_p_values_icr_enabled)
colnames(ordered_p_values_icr_enabled) <- c("MR","Adj_Pval")
ordered_p_values_icr_enabled$MR <- as.character(as.vector(ordered_p_values_icr_enabled$MR))
ordered_p_values_icr_enabled$Adj.Pval <- as.numeric(as.vector(ordered_p_values_icr_enabled$Adj_Pval))
write.table(ordered_p_values_icr_enabled,"../Results/Revised_Text_Results/Ordered_Pvalues_ICR_Enabled_MRs.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get activity of MRs specific to ICR Low of interest from enabled cancers
mrs_of_interest_enabled <- icr_enabled_median_comparison[icr_enabled_median_comparison$FC_Median<0 &
                                                         icr_enabled_median_comparison$Pval<0.05,]$MR
enabled_activity_df <- NULL
for (cancer in icr_enabled)
{
  load(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  
  high_indices_table <- read.table(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_low_indices.csv"),header=TRUE)
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
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ TopMR, nrow=5, ncol=7) + xlab("Cancer") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top MRs of specific to ICR Low in ICR Enabled Cancers (>=4 out of 8)") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_4.jpg",plot = supp_p3, device = jpeg(), width = 14, height=12, units = "in", dpi = 300)
dev.off()

icr_enabled_median_comparison <- icr_enabled_median_comparison[order(icr_enabled_median_comparison$FC_Median,decreasing = T),]
final_icr_enabled_median_comparison <- icr_enabled_median_comparison[icr_enabled_median_comparison$Pval<0.05,]
write.table(final_icr_enabled_median_comparison,file="../Results/Revised_Text_Results/All_MRs_TFs_ICR_Enabled_Supplementary_Table_Revised(v2)_S6b.csv",
            row.names = F, col.names = T, sep = " & ", quote=F)

#Perform Analysis for ICR Disabled
###########################################################################################################################
tfs_disabled_activity_info <- get_tf_median_activity(icr_disabled[1])
for (i in 2:length(icr_disabled))
{
  cancer <- icr_disabled[i]
  temp <- get_tf_median_activity(cancer)
  tfs_disabled_activity_info <- rbind(tfs_disabled_activity_info,temp)
}
tfs_disabled_activity_info <- as.data.frame(tfs_disabled_activity_info)
colnames(tfs_disabled_activity_info) <- c("TFs","Cancer","Median_ICR_High","Median_ICR_Low")
tfs_disabled_activity_info$TFs <- as.character(as.vector(tfs_disabled_activity_info$TFs))
tfs_disabled_activity_info$Cancer <- as.character(as.vector(tfs_disabled_activity_info$Cancer))
tfs_disabled_activity_info$Median_ICR_High <- as.numeric(as.vector(tfs_disabled_activity_info$Median_ICR_High))
tfs_disabled_activity_info$Median_ICR_Low <- as.numeric(as.vector(tfs_disabled_activity_info$Median_ICR_Low))

#Select those transcription regulators which are TFs with regulons of size >= 10 in all the ICR Disabled Cancers
common_tfs_disabled_list <- names(which(table(tfs_disabled_activity_info$TFs)>3))
common_tfs_disabled_activity_info <- tfs_disabled_activity_info[tfs_disabled_activity_info$TFs %in% common_tfs_disabled_list,]

first_cancer <- icr_disabled[1];
list_common_mrs_disabled <- get_common_mrs(first_cancer)
topMRs_disabled <- rep(list("LGG"),length(list_common_mrs_disabled))
names(topMRs_disabled) <- list_common_mrs_disabled
for(i in 2:length(icr_disabled))
{
  cancer_type <- icr_disabled[i]
  list_common_mrs_disabled <- get_common_mrs(cancer_type)
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

temp_disabled_df <- disabled_df[disabled_df$No_ICR_Cancers>=2,]
temp_disabled_df <- temp_disabled_df[temp_disabled_df$MR %in% common_tfs_disabled_list, ]
write.table(temp_disabled_df,"../Results/Revised_Text_Results/All_MRs_TFs_ICR_Disabled_Supplementary_Table_Revised(v2)_S7.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get list of common MRs which are also TFs in 2 out of 4 cancers
############################################################################################
common_mrs_disabled_list <- temp_disabled_df$MR;
common_mrs_disabled_list <- sort(common_mrs_disabled_list)
common_mrs_disabled_activity_matrix <- matrix(0,nrow=length(common_mrs_disabled_list),ncol=2*length(icr_disabled))
rownames(common_mrs_disabled_activity_matrix) <- common_mrs_disabled_list
colnames(common_mrs_disabled_activity_matrix) <- c(icr_disabled,icr_disabled)
for (i in 1:length(icr_disabled))
{
  cancer <- icr_disabled[i]
  temp_df <- common_tfs_disabled_activity_info[common_tfs_disabled_activity_info$TFs %in% common_mrs_disabled_list 
                                              & common_tfs_disabled_activity_info$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$TFs),]
  icr_high <- temp_df$Median_ICR_High
  icr_low <- temp_df$Median_ICR_Low
  common_mrs_disabled_activity_matrix[common_mrs_disabled_list,i] <- icr_high
  common_mrs_disabled_activity_matrix[common_mrs_disabled_list,(i+length(icr_disabled))] <- icr_low
}

colcol <- matrix(0,nrow=2*length(icr_disabled),ncol=1)
colcol[c(1:length(icr_disabled)),1] <- "yellow"
colcol[c((length(icr_disabled)+1):(2*length(icr_disabled))),1] <- "green"

#Perform Wilcox ranksum test or Mann-Whitney test to identify the MRs whose activity between the ICR High and ICR Low are significant for ICR Disabled cancers
#######################################################################################################
output_df <- common_tfs_disabled_activity_info
nes_info <- NULL
nes_pval <- NULL
for (i in 1:length(common_mrs_disabled_list))
{
  mr <- common_mrs_disabled_list[i]
  mr_icr_high_values <- NULL
  mr_icr_low_values <- NULL
  for (cancer in icr_disabled)
  {
    mr_activity_list <- get_tf_activity_high_low(cancer,mr)
    mr_icr_high_values <- c(mr_icr_high_values,mr_activity_list[[1]])
    mr_icr_low_values <- c(mr_icr_low_values,mr_activity_list[[2]])
  }
  nes_info <- c(nes_info,median(mr_icr_high_values)-median(mr_icr_low_values))
  nes_pval <- c(nes_pval,wilcox.test(mr_icr_high_values,mr_icr_low_values,exact=F)$p.value)
  
}
icr_disabled_median_comparison <- cbind(common_mrs_disabled_list,nes_info,nes_pval)
icr_disabled_median_comparison <- as.data.frame(icr_disabled_median_comparison)
colnames(icr_disabled_median_comparison) <- c("MR","FC_Median","Pval")
icr_disabled_median_comparison$MR <- as.character(as.vector(icr_disabled_median_comparison$MR))
icr_disabled_median_comparison$FC_Median <- round(as.numeric(as.vector(icr_disabled_median_comparison$FC_Median)),3)
icr_disabled_median_comparison$Pval <- p.adjust(as.numeric(as.vector(icr_disabled_median_comparison$Pval)),method = "fdr")
final_common_mrs_disabled_list <- icr_disabled_median_comparison[icr_disabled_median_comparison$Pval<0.05,]$MR

final_common_mrs_disabled_activity_matrix <- common_mrs_disabled_activity_matrix[final_common_mrs_disabled_list,]

#Figure 4B
#=========================================================================================================
pdf("../Results/Revised_Figures/Pdfs/Common_TopMRs_Activity_ICR_Disabled_Figure_Revised(v2)_4B.pdf",height = 12, width=14, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.65)
p2 <- heatmap.3(final_common_mrs_disabled_activity_matrix, Rowv = TRUE, Colv=TRUE, col = bluered(100), scale="none", main= "Median Activity of Common MRs in ICR Disabled Cancers", # (>=2 out of 4)",
                dendrogram = "both", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.6, cexCol = 2, cellnote = ifelse(final_common_mrs_disabled_activity_matrix>0,"+","-"), notecex = 0.75, notecol = "black")
dev.off()

ordered_mrs_disabled <- rev(rownames(final_common_mrs_disabled_activity_matrix)[p2$rowInd])
ordered_p_values_icr_disabled <- NULL
for (i in 1:length(ordered_mrs_disabled))
{
  mr <- ordered_mrs_disabled[i]
  adj_pval <- icr_disabled_median_comparison[icr_disabled_median_comparison$MR==mr,]$Pval
  temp <- cbind(mr,adj_pval)
  ordered_p_values_icr_disabled <- rbind(ordered_p_values_icr_disabled,temp)
}
ordered_p_values_icr_disabled <- as.data.frame(ordered_p_values_icr_disabled)
colnames(ordered_p_values_icr_disabled) <- c("MR","Adj_Pval")
ordered_p_values_icr_disabled$MR <- as.character(as.vector(ordered_p_values_icr_disabled$MR))
ordered_p_values_icr_disabled$Adj.Pval <- as.numeric(as.vector(ordered_p_values_icr_disabled$Adj_Pval))
write.table(ordered_p_values_icr_disabled,"../Results/Revised_Text_Results/Ordered_Pvalues_ICR_Disabled_MRs.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get activity of mrs of interest from ICR Disabled cancers
mrs_of_interest_disabled <- icr_disabled_median_comparison[icr_disabled_median_comparison$FC_Median<0 &
                                                           icr_disabled_median_comparison$Pval<0.05,]$MR
disabled_activity_df <- NULL
for (cancer in icr_disabled)
{
  load(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  
  high_indices_table <- read.table(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_low_indices.csv"),header=TRUE)
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
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ TopMR, nrow=10, ncol=9) + xlab("Cancer") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top MRs of specific to ICR Low in ICR Enabled Cancers (>=2 out of 4)") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_5.jpg",plot = supp_p4, device = jpeg(), width = 14, height=12, units = "in", dpi = 300)
dev.off()

icr_disabled_median_comparison <- icr_disabled_median_comparison[order(icr_disabled_median_comparison$FC_Median,decreasing = T),]
final_icr_disabled_median_comparison <- icr_disabled_median_comparison[icr_disabled_median_comparison$Pval<0.05,]
write.table(final_icr_disabled_median_comparison,file="../Results/Revised_Text_Results/All_MRs_TFs_ICR_Disabled_Supplementary_Table_Revised(v2)_S7b.csv",
            row.names = F, col.names = T, sep = " & ", quote=F)
#=========================================================================================================================

#Activity of MRs specific to ICR High in all 12 cancers of interest (Figures 4C, 4E)
###########################################################################################################################
mrs_icr_high_all_cancers <- union(final_icr_enabled_median_comparison[final_icr_enabled_median_comparison$FC_Median>0,]$MR,
                                  final_icr_disabled_median_comparison[final_icr_disabled_median_comparison$FC_Median>0,]$MR)

all_cancers <- c(icr_enabled,icr_disabled)
mrs_icr_high_activity_matrix <- matrix(0,nrow=length(mrs_icr_high_all_cancers),ncol=(length(icr_enabled)+length(icr_disabled)))
rownames(mrs_icr_high_activity_matrix) <- mrs_icr_high_all_cancers
colnames(mrs_icr_high_activity_matrix) <- all_cancers
mrs_high_df <- NULL
for (i in 1:length(mrs_icr_high_all_cancers))
{
  topmr <- mrs_icr_high_all_cancers[i]
  icr_enabled_high_activity_values <- NULL
  icr_disabled_high_activity_values <- NULL
  for (j in 1:length(all_cancers))
  {
    cancer <- all_cancers[j]
    if (cancer %in% icr_enabled)
    {
      out_activity <- get_tf_activity_high_low(cancer,topmr)
      icr_enabled_high_activity_values <- c(icr_enabled_high_activity_values,out_activity[[1]])
      mrs_icr_high_activity_matrix[topmr,cancer] <- median(out_activity[[1]])
    } 
    else if (cancer %in% icr_disabled)
    {
      out_activity <- get_tf_activity_high_low(cancer,topmr)
      icr_disabled_high_activity_values <- c(icr_disabled_high_activity_values,out_activity[[1]])
      mrs_icr_high_activity_matrix[topmr,cancer] <- median(out_activity[[1]])
    }
  }
  median_icr_enabled_high <- median(icr_enabled_high_activity_values)
  median_icr_disabled_high <- median(icr_disabled_high_activity_values)
  fc_info <- median(icr_enabled_high_activity_values)-median(icr_disabled_high_activity_values)
  pval <- wilcox.test(icr_enabled_high_activity_values,icr_disabled_high_activity_values,exact=F)$p.value
  temp <- cbind(topmr,fc_info,median_icr_enabled_high,median_icr_disabled_high,pval)
  mrs_high_df <- rbind(mrs_high_df,temp)
}
mrs_high_df <- as.data.frame(mrs_high_df)
colnames(mrs_high_df) <- c("MR","FC_Median","Median_ICR_Enabled","Median_ICR_Disabled","Padj")
mrs_high_df$MR <- as.character(as.vector(mrs_high_df$MR))
mrs_high_df$Padj <- p.adjust(as.numeric(as.vector(mrs_high_df$Padj)),method="fdr")
mrs_high_df$FC_Median <- round(as.numeric(as.vector(mrs_high_df$FC_Median)),3)
mrs_high_df$Median_ICR_Enabled <- round(as.numeric(as.vector(mrs_high_df$Median_ICR_Enabled)),3)
mrs_high_df$Median_ICR_Disabled <- round(as.numeric(as.vector(mrs_high_df$Median_ICR_Disabled)),3)

top_positive_common_mrs_high <- mrs_high_df[mrs_high_df$Median_ICR_Enabled>=0 & mrs_high_df$Median_ICR_Disabled>=0,]$MR
top_negative_enabled_positive_disabled_mrs_high <- mrs_high_df[mrs_high_df$Median_ICR_Enabled<0.0 & 
                                                               mrs_high_df$Median_ICR_Disabled>0 &
                                                               mrs_high_df$Padj<0.05,]$MR


interesting_mr_high <- c(top_positive_common_mrs_high,top_negative_enabled_positive_disabled_mrs_high)
interesting_mr_high_df <- mrs_high_df[mrs_high_df$MR %in% interesting_mr_high,]
interesting_mr_high_df <- interesting_mr_high_df[order(interesting_mr_high_df$Median_ICR_Enabled,decreasing=T),]
write.table(interesting_mr_high_df,"../Results/Revised_Text_Results/All_MRS_ICR_High_Supplementary_Table_S7c.csv",row.names=F,col.names=T,sep="&", quote=F)


positive_mrs_icr_high_activity_matrix <- mrs_icr_high_activity_matrix[top_positive_common_mrs_high,]
pdne_mrs_icr_high_activity_matrix <- mrs_icr_high_activity_matrix[top_negative_enabled_positive_disabled_mrs_high,] 

colcol <- matrix(0,nrow=length(all_cancers),ncol=1)
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_cancers)),1] <- "#660066"

#Figure 4C
#====================================================================================
pdf("../Results/Revised_Figures/Pdfs/Common_TopMRs_Activity_ICR_High_Figure_Revised(v2)_4C.pdf",height = 10, width=14, pointsize = 12)
p3 <- heatmap.3(positive_mrs_icr_high_activity_matrix, Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Median Activity of MRs specific to ICR High Phenotype for 12 cancers in ICR High samples",
                dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.5, cexCol = 1.5, cellnote = ifelse(positive_mrs_icr_high_activity_matrix>=0,"+","-"), notecex = 1.0, notecol = "black")
dev.off()

hc.rows <- hclust(dist(pdne_mrs_icr_high_activity_matrix),method='ward.D2')
colcol <- matrix(0,nrow=length(all_cancers),ncol=2)
colnames(colcol) <- c("Enabled/Disabled","ICR High")
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_cancers)),1] <- "#660066"
colcol[c(1:12),2] <- "yellow"

#Figure 4E
#====================================================================================
pdf("../Results/Revised_Figures/Pdfs/Different_TopMRs_Activity_ICR_High_Figure_Revised(v2)_4E.pdf",height = 10, width=12, pointsize = 12)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.0)
heatmap.3(pdne_mrs_icr_high_activity_matrix[hc.rows$order,], Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Median Activity of ICR High MRs different between ICR Enabled & ICR Disabled cancers",
          dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 3,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.5, cexCol = 1.5, cellnote = ifelse(pdne_mrs_icr_high_activity_matrix[hc.rows$order,]<=0.0,"-","+"), notecex = 2, notecol = "black")
dev.off()
#==================================================================================================================================

#Activity of MRs specific to ICR High in all 12 cancers of interest (Figures 4D, 4F)
#==============================================================================================================================
mrs_icr_low_all_cancers <- union(final_icr_enabled_median_comparison[final_icr_enabled_median_comparison$FC_Median<0,]$MR,
                                  final_icr_disabled_median_comparison[final_icr_disabled_median_comparison$FC_Median<0,]$MR)

all_cancers <- c(icr_enabled,icr_disabled)
mrs_icr_low_activity_matrix <- matrix(0,nrow=length(mrs_icr_low_all_cancers),ncol=(length(icr_enabled)+length(icr_disabled)))
rownames(mrs_icr_low_activity_matrix) <- mrs_icr_low_all_cancers
colnames(mrs_icr_low_activity_matrix) <- all_cancers
mrs_low_df <- NULL
for (i in 1:length(mrs_icr_low_all_cancers))
{
  topmr <- mrs_icr_low_all_cancers[i]
  icr_enabled_low_activity_values <- NULL
  icr_disabled_low_activity_values <- NULL
  for (j in 1:length(all_cancers))
  {
    cancer <- all_cancers[j]
    if (cancer %in% icr_enabled)
    {
      out_activity <- get_tf_activity_high_low(cancer,topmr)
      icr_enabled_low_activity_values <- c(icr_enabled_low_activity_values,out_activity[[2]])
      mrs_icr_low_activity_matrix[topmr,cancer] <- median(out_activity[[2]])
    } 
    else if (cancer %in% icr_disabled)
    {
      out_activity <- get_tf_activity_high_low(cancer,topmr)
      icr_disabled_low_activity_values <- c(icr_disabled_low_activity_values,out_activity[[2]])
      mrs_icr_low_activity_matrix[topmr,cancer] <- median(out_activity[[2]])
    }
  }
  median_icr_enabled_low <- median(icr_enabled_low_activity_values)
  median_icr_disabled_low <- median(icr_disabled_low_activity_values)
  fc_info <- median(icr_enabled_low_activity_values)-median(icr_disabled_low_activity_values)
  pval <- wilcox.test(icr_enabled_low_activity_values,icr_disabled_low_activity_values,exact=F)$p.value
  temp <- cbind(topmr,fc_info,median_icr_enabled_low,median_icr_disabled_low,pval)
  mrs_low_df <- rbind(mrs_low_df,temp)
}
mrs_low_df <- as.data.frame(mrs_low_df)
colnames(mrs_low_df) <- c("MR","FC_Median","Median_ICR_Enabled","Median_ICR_Disabled","Padj")
mrs_low_df$MR <- as.character(as.vector(mrs_low_df$MR))
mrs_low_df$Padj <- p.adjust(as.numeric(as.vector(mrs_low_df$Padj)),method="fdr")
mrs_low_df$FC_Median <- round(as.numeric(as.vector(mrs_low_df$FC_Median)),3)
mrs_low_df$Median_ICR_Enabled <- round(as.numeric(as.vector(mrs_low_df$Median_ICR_Enabled)),3)
mrs_low_df$Median_ICR_Disabled <- round(as.numeric(as.vector(mrs_low_df$Median_ICR_Disabled)),3)

top_positive_common_mrs_low <- mrs_low_df[mrs_low_df$Median_ICR_Enabled>0 & mrs_low_df$Median_ICR_Disabled>0,]$MR
top_positive_enabled_negative_disabled_mrs_low <- mrs_low_df[mrs_low_df$Median_ICR_Enabled>=0.0 & 
                                                                 mrs_low_df$Median_ICR_Disabled<0 &
                                                                 mrs_low_df$Padj<0.05,]$MR

interesting_mr_low <- c(top_positive_common_mrs_low,top_positive_enabled_negative_disabled_mrs_low)
interesting_mr_low_df <- mrs_low_df[mrs_low_df$MR %in% interesting_mr_low,]
interesting_mr_low_df <- interesting_mr_low_df[order(interesting_mr_low_df$Median_ICR_Disabled,decreasing=T),]
write.table(interesting_mr_low_df,"../Results/Revised_Text_Results/All_MRS_ICR_Low_Supplementary_Table_S7d.csv",row.names=F,col.names=T,sep="&", quote=F)

positive_mrs_icr_low_activity_matrix <- mrs_icr_low_activity_matrix[top_positive_common_mrs_low,]
nepd_mrs_icr_low_activity_matrix <- mrs_icr_low_activity_matrix[top_positive_enabled_negative_disabled_mrs_low,]

colcol <- matrix(0,nrow=length(all_cancers),ncol=1)
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_cancers)),1] <- "#660066"

#Figure 4D
#====================================================================================
pdf("../Results/Revised_Figures/Pdfs/Common_TopMRs_Activity_ICR_Low_Figure_Revised(v2)_4D.pdf",height = 10, width=14, pointsize = 12)
p3 <- heatmap.3(positive_mrs_icr_low_activity_matrix, Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Median Activity of MRs specific to ICR Low Phenotype for 12 cancers in ICR High samples",
                dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.65, cexCol = 1.5, cellnote = ifelse(positive_mrs_icr_low_activity_matrix>=0,"+","-"), notecex = 1.0, notecol = "black")
dev.off()

#Plot 4F
hc.rows <- hclust(dist(nepd_mrs_icr_low_activity_matrix),method='ward.D2')
colcol <- matrix(0,nrow=length(all_cancers),ncol=2)
colnames(colcol) <- c("Enabled/Disabled","ICR High")
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_cancers)),1] <- "#660066"
colcol[c(1:12),2] <- "green"

pdf("../Results/Revised_Figures/Pdfs/Different_TopMRs_Activity_ICR_Low_Figure_Revised(v2)_4F.pdf",height = 10, width=12, pointsize = 12)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.0)
heatmap.3(nepd_mrs_icr_low_activity_matrix[hc.rows$order,], Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Median Activity of ICR Low MRs different between ICR Enabled & ICR Disabled cancers",
          dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 3,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.5, cexCol = 1.5, cellnote = ifelse(nepd_mrs_icr_low_activity_matrix[hc.rows$order,]<0.0,"-","+"), notecex = 2, notecol = "black")
dev.off()

#Make figures for downstream enrichment analysis (Supp 6A, 6B)
##################################################################################################################################
icr_high_goterms_df <- read.table("../Results/Revised_Text_Results/Final_Enriched_GOTerms_ICR_High.csv",header=TRUE,sep="\t")
temp_df <- icr_high_goterms_df[-log10(icr_high_goterms_df$p.value)>40,]

supp_p6a <- ggplot(icr_high_goterms_df,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) + geom_point(aes(color=as.factor(term_level))) + 
  xlab("Percentage of MRs in each GO Term") + ylab("-log10(Pvalue)") + scale_shape_manual(name = "GO Terms", values = 1:nlevels(as.factor(icr_high_goterms_df$term_category))) +
  scale_color_manual(name = "Term Level", values =  c("grey","violet","orange")) +
  geom_hline(yintercept=c(0,20,40), color=c("black","blue","red")) + geom_vline(xintercept = 0) + 
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlim(c(0.0,80))+
  geom_text(aes(x=100*generatio,y=-log10(p.value),label=term_name),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.025, angle = 0,
            size = 3, check_overlap = T) +
  ggtitle("Enriched GO Terms for ICR High Phenotype") + theme(text = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_6A.jpg",plot = supp_p6a, device = jpeg(), width = 7, height=6, units = "in", dpi = 300)
dev.off()

icr_low_goterms_df <- read.table("../Results/Revised_Text_Results/Final_Enriched_GOTerms_ICR_Low.csv",header=TRUE,sep="\t")
temp_df <- icr_low_goterms_df[-log10(icr_low_goterms_df$p.value)>20,]

supp_p6b <- ggplot(icr_low_goterms_df,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) + geom_point(aes(color=as.factor(term_level))) + 
  xlab("Percentage of MRs in each GO Term") + ylab("-log10(Pvalue)") + scale_shape_manual(name = "GO Terms", values = 1:nlevels(as.factor(icr_low_goterms_df$term_category))) +
  scale_color_manual(name = "Term Level", values =  c("grey","violet","orange")) +
  geom_hline(yintercept=c(0,10,20), color=c("black","blue","red")) + geom_vline(xintercept = 0) + 
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(aes(x=100*generatio,y=-log10(p.value),label=term_name),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.025, angle = 0,
            size = 3, check_overlap = T) +
  coord_cartesian(xlim =c(0, 15), ylim = c(0, 30)) +
  ggtitle("Enriched GO Terms for ICR Low Phenotype") + theme(text = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_6B.jpg",plot = supp_p6b, device = jpeg(), width = 7, height=6, units = "in", dpi = 300)
dev.off()

#Make Figure 5A_2
#===========================================================================================================
colors <- c('#3cb44b', '#4363d8', '#f58231', '#911eb4', '#f032e6', '#46f0f0' , '#bcf60c', '#a9a9a9', '#e6194b', '#ffe119', '#ffe119', '#e6beff')
icr_low_pathways_df <- read.table("../Results/Revised_Text_Results/Final_Enriched_Pathways_ICR_Low_v2.csv",header=TRUE,sep="\t")
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
ggsave(filename = "../Results/Revised_Figures/Pdfs/Enriched_Pathways_ICR_Low_Figure_Revised(v2)_5A_2.pdf",plot=p5,device=pdf(), height=12, width=15, units="in", dpi=300)
dev.off()

#Supplementary Figure S7.A
#================================================================================================================
icr_high_pathways_df <- read.table("../Results/Revised_Text_Results/Final_Enriched_Pathways_ICR_High_v2.csv",header=TRUE,sep="\t")
icr_high_pathways_df$Description <- as.character(as.vector(icr_high_pathways_df$Description))
icr_high_pathways_df$GeneRatio <- as.character(as.vector(icr_high_pathways_df$GeneRatio))
icr_high_pathways_df$Description <- paste0(icr_high_pathways_df$Description," [",icr_high_pathways_df$GeneRatio,"] ")
icr_high_pathways_df$genes <- as.character(as.vector(icr_high_pathways_df$genes))
icr_high_pathways_df <- icr_high_pathways_df[order(icr_high_pathways_df$generatio),]
supp_p7a <- ggplot(data=icr_high_pathways_df,aes(x=generatio,y=reorder(Description,generatio),size=-log10(pvalue)))+
  geom_point(aes(color=-log10(pvalue)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(name = "-log10(P.adjust)")+ guides(color=guide_legend(), size = guide_legend())+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Enriched Pathways for ICR High Phenotype") + theme(text = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y=element_text(color=colors[as.factor(icr_high_pathways_df$clusters)]))
ggsave(filename = "../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_7A.jpg",plot=supp_p7a,device=jpeg(), height=15, width=12, units="in", dpi=300)
dev.off()

#Supplementary Figure S7B
#=====================================================================================================================
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.1), c("blue", "white", "red"))
col_fun(seq(-3, 3))

icr_high_pathways_df <- icr_high_pathways_df[order(icr_high_pathways_df$generatio,decreasing = T),]
icr_high_pathways_gene_involved <- matrix(0,nrow=nrow(icr_high_pathways_df),ncol=length(rownames(positive_mrs_icr_high_activity_matrix)))
rownames(icr_high_pathways_gene_involved) <- icr_high_pathways_df$Description
colnames(icr_high_pathways_gene_involved) <- rownames(positive_mrs_icr_high_activity_matrix)

for (i in 1:nrow(icr_high_pathways_df))
{
  genes_involved <- unlist(strsplit(icr_high_pathways_df[i,]$genes,split="; "))
  median_activity_scores <- rowMedians(positive_mrs_icr_high_activity_matrix[genes_involved,])
  icr_high_pathways_gene_involved[i,genes_involved] <- median_activity_scores
}

jpeg("../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_7B.jpg",height = 900, width=1500, units="px",pointsize = 12)
Heatmap(icr_high_pathways_gene_involved,cluster_rows = FALSE, cluster_columns = FALSE, col=col_fun, row_names_side = "left", 
        name="Median Activity", row_names_gp = gpar(fontsize = 6, col=colors[as.factor(icr_high_pathways_df$clusters)]), column_names_gp = gpar(fontsize = 8), 
        column_title = "Median Activity of MRs in Enriched Pathways specific to ICR High")
dev.off()

#Figure 5A_1
#================================================================================================================
rev_icr_low_pathways_df <- icr_low_pathways_df[order(icr_low_pathways_df$generatio,decreasing = T),]
icr_low_pathways_gene_involved <- matrix(0,nrow=nrow(rev_icr_low_pathways_df),ncol=length(rownames(positive_mrs_icr_low_activity_matrix)))
rownames(icr_low_pathways_gene_involved) <- rev_icr_low_pathways_df$Description
colnames(icr_low_pathways_gene_involved) <- rownames(positive_mrs_icr_low_activity_matrix)
for (i in 1:nrow(rev_icr_low_pathways_df))
{
  genes_involved <- unlist(strsplit(rev_icr_low_pathways_df[i,]$genes,split="; "))
  median_activity_scores <- rowMedians(positive_mrs_icr_low_activity_matrix[genes_involved,])
  icr_low_pathways_gene_involved[i,genes_involved] <- median_activity_scores
}

#Make the Sankey plot
icr_low_pathways_gene_involved_to_consider <- icr_low_pathways_gene_involved[,which(colSums(icr_low_pathways_gene_involved)>0)]

#Figure S7C
jpeg("../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_7C.jpg",height = 600, width=900, units="px",pointsize = 15, res = 0.75)
op <- par(oma=c(10,7,1,1))
Heatmap(icr_low_pathways_gene_involved_to_consider,cluster_rows = FALSE, cluster_columns = FALSE, col=col_fun, row_names_side = "left", 
        name="Activity", row_names_gp = gpar(fontsize = 6, col=colors[as.factor(rev_icr_low_pathways_df$clusters)]), column_names_gp = gpar(fontsize = 10), 
        column_title = "Median Activity of MRs in Enriched Pathways specific to ICR Low")
dev.off()

ICR_Low <- list()
ICR_Low$nodes <- data.frame(name = c(colnames(icr_low_pathways_gene_involved_to_consider),rownames(icr_low_pathways_gene_involved_to_consider)))
nodesgroup <- c(rep("white",length(colnames(icr_low_pathways_gene_involved_to_consider))),colors[rev_icr_low_pathways_df$clusters])
ICR_Low$nodes$nodesgroup <- nodesgroup
edgelist <- NULL
for (i in 1:nrow(icr_low_pathways_gene_involved_to_consider))
{
  for (j in 1:ncol(icr_low_pathways_gene_involved_to_consider))
  {
    pathway <- rownames(icr_low_pathways_gene_involved_to_consider)[i]
    mr <- colnames(icr_low_pathways_gene_involved_to_consider)[j]
    if (icr_low_pathways_gene_involved_to_consider[pathway,mr]>0)
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
  range = c(rep("white",length(colnames(icr_low_pathways_gene_involved_to_consider))),colors[rev_icr_low_pathways_df$clusters]),
  domain = ICR_Low$nodes$name,
  nodes = ICR_Low$nodes,
  stringsAsFactors = FALSE
)

p <- sankeyNetwork(Links = ICR_Low$links, Nodes = ICR_Low$nodes, Source = "source",
                   Target = "target", Value = "value", LinkGroup = "type", NodeID = "name", NodeGroup = "nodesgroup",
                   units = "", fontSize = 16, nodeWidth = 25, iterations=0, fontFamily = "Arial", 
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);")) 
#widget2png(p, "../Results/Revised_Figures/Svgs_Jpgs/Sankey_Plot_ICR_Low_Figure_Revised(v2)_5A.png")


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
#     plot_bgcolor = 'white',
#     paper_bgcolor = 'white'
#   )

# plotly_IMAGE(p, format = "png", file = "Results/Figures/Pdfs/Sankey_ICR_Low.png", width = 600, height = 480)


# #Load the images 4A and 4B
im1 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Common_TopMRs_Activity_ICR_Enabled_Figure_Revised(v2)_4A.svg")
p1 <- pictureGrob(im1)
 
im2 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Common_TopMRs_Activity_ICR_Disabled_Figure_Revised(v2)_4B.svg")
p2 <- pictureGrob(im2)
 
im3 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Common_TopMRs_Activity_ICR_High_Figure_Revised(v2)_4C.svg")
p3 <- pictureGrob(im3)
 
im4 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Common_TopMRs_Activity_ICR_Low_Figure_Revised(v2)_4D.svg")
p4 <- pictureGrob(im4)
 
im5 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Different_TopMRs_Activity_ICR_High_Figure_Revised(v2)_4E.svg")
p5 <- pictureGrob(im5)
 
im6 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Different_TopMRs_Activity_ICR_Low_Figure_Revised(v2)_4F.svg")
p6 <- pictureGrob(im6)
 
g_final <- ggarrange(p1,p2,p3,p4,p5,p6,nrow=2,ncol=3,labels=c("A","B","C","D","E","F"),widths=c(2,2))
ggsave(filename = paste0("../Results/Revised_Figures/Figure_Revised(v2)_4.pdf"),device = pdf(), plot = g_final, dpi = 300, width = 16, height= 14, units = "in" )
dev.off()
