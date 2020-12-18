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
library(gridExtra)
library(ggpubr)
library(grImport2)
library(matrixStats)
warnings("off")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

registerDoMC(20)

setwd('/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Old_Files/ICR_All_Info/scripts/')
setwd('.')

source('gene-reverse-network.R')
source('get_functions.R')

revise_wilcox_test_info <- function(wilcox_test_info)
{
  wilcox_test_info$Pval <- round(wilcox_test_info$Pval,5)
  wilcox_test_info$Mean1 <- round(wilcox_test_info$Mean1,3)
  wilcox_test_info$Mean2 <- round(wilcox_test_info$Mean2,3)
  wilcox_test_info$FC_Mean <- round(wilcox_test_info$FC_Mean,3)
  return(wilcox_test_info)
}

icr_enabled <- c("BLCA","BRCA","HNSC","LIHC","SARC","SKCM","STAD","UCEC")
icr_disabled <- c("LGG","KIRC","PAAD","UVM")

common_mrs_enabled_list <- NULL
common_mrs_disabled_list <- NULL

output_df <- read.table("../Results/Revised_Text_Results/All_TF_Activity_Information_v2.csv",header = TRUE)
output_df$MR <- as.character(as.vector(output_df$MR))
output_df$Cancer <- as.character(as.vector(output_df$Cancer))
colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#a9a9a9', '#008080', '#e6beff')

common_mrs_enabled_list <- c(output_df[output_df$Cancer==icr_enabled[1],]$MR)
for (i in 2:length(icr_enabled))
{
  new_mrs <- output_df[output_df$Cancer==icr_enabled[i],]$MR
  common_mrs_enabled_list <- intersect(common_mrs_enabled_list,new_mrs)
}

common_mrs_enabled_activity_matrix <- matrix(0,nrow=length(common_mrs_enabled_list),ncol=2*length(icr_enabled))
common_mrs_enabled_list <- sort(common_mrs_enabled_list)
rownames(common_mrs_enabled_activity_matrix) <- common_mrs_enabled_list
colnames(common_mrs_enabled_activity_matrix) <- c(icr_enabled,icr_enabled)
for (i in 1:length(icr_enabled))
{
  cancer <- icr_enabled[i]
  temp_df <- output_df[output_df$MR %in% common_mrs_enabled_list & output_df$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$MR),]
  icr_high <- temp_df$Median_ICR_High
  icr_low <- temp_df$Median_ICR_Low
  common_mrs_enabled_activity_matrix[common_mrs_enabled_list,i] <- icr_high
  common_mrs_enabled_activity_matrix[common_mrs_enabled_list,(i+length(icr_enabled))] <- icr_low
}

colcol <- matrix(0,nrow=2*length(icr_enabled),ncol=1)
colcol[c(1:length(icr_enabled)),1] <- "yellow"
colcol[c((length(icr_enabled)+1):(2*length(icr_enabled))),1] <- "green"

#Figure 3A 
#=========================================================================================================
pdf("../Results/Revised_Figures/Pdfs/Common_TopMRs_Activity_ICR_Enabled_Figure_Revised(v2)_3A.pdf",height = 10, width=12, pointsize = 14)
heatmap.3(common_mrs_enabled_activity_matrix, Rowv = TRUE, Colv=TRUE, col = bluered(100), scale="none", main= "Median Activity of MRs common across all ICR Enabled Cancers",
          labRow = rownames(common_mrs_enabled_activity_matrix), labCol = colnames(common_mrs_enabled_activity_matrix), dendrogram = "col", 
          key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.25, cexCol = 2, cellnote = ifelse(common_mrs_enabled_activity_matrix>0,"+","-"), notecex = 2, notecol = "black")
dev.off()

#Supplementary information about top common MRs for ICR ENabled cancer based on fold change values
A <- common_mrs_enabled_activity_matrix[common_mrs_enabled_list,c(1:length(icr_enabled))];
B <- common_mrs_enabled_activity_matrix[common_mrs_enabled_list,c((length(icr_enabled)+1):(2*length(icr_enabled)))]
wilcox_test_info <- perform_wilcox_test(A,B)
wilcox_test_info <- wilcox_test_info[order(wilcox_test_info$FC_Mean,decreasing=TRUE),]
wilcox_test_info <- revise_wilcox_test_info(wilcox_test_info)
write.table(wilcox_test_info,"../Results/Revised_Text_Results/Supplementary_Table_S4_Common_MRs_ICR_Enabled_Revised(v2).csv",row.names=T,col.names=T,quote=F)
top_10_common_mrs_enabled <- rownames(wilcox_test_info)[1:10]

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
  
  for (topmr in top_10_common_mrs_enabled)
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

p3 <- ggplot(data = enabled_activity_df, aes(x=TopMR, y=Activity_Value)) + 
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ Cancer, nrow=2, ncol=4) + xlab("Top Common MRs") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top 10 MRs common to ICR Enabled Cancers ") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))

ggsave(file="../Results/Revised_Figures/Pdfs/Common_TopMRs_Activity_Figure_Revised(v2)_3C.pdf",plot = p3, device = pdf(), width = 12, height=12, units = "in", dpi = 300)
dev.off()

#Perform the same analysis for ICR Disabled cancers
#==============================================================================================================================
common_mrs_disabled_list <- c(output_df[output_df$Cancer==icr_disabled[1],]$MR)
for (i in 2:length(icr_disabled))
{
  new_mrs <- output_df[output_df$Cancer==icr_disabled[i],]$MR
  common_mrs_disabled_list <- intersect(common_mrs_disabled_list,new_mrs)
}

common_mrs_disabled_activity_matrix <- matrix(0,nrow=length(common_mrs_disabled_list),ncol=2*length(icr_disabled))
common_mrs_disabled_list <- sort(common_mrs_disabled_list)
rownames(common_mrs_disabled_activity_matrix) <- common_mrs_disabled_list
colnames(common_mrs_disabled_activity_matrix) <- c(icr_disabled,icr_disabled)
for (i in 1:length(icr_disabled))
{
  cancer <- icr_disabled[i]
  temp_df <- output_df[output_df$MR %in% common_mrs_disabled_list & output_df$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$MR),]
  icr_high <- temp_df$Median_ICR_High
  icr_low <- temp_df$Median_ICR_Low
  common_mrs_disabled_activity_matrix[common_mrs_disabled_list,i] <- icr_high
  common_mrs_disabled_activity_matrix[common_mrs_disabled_list,(i+length(icr_disabled))] <- icr_low
}

colcol <- matrix(0,nrow=2*length(icr_disabled),ncol=1)
colcol[c(1:length(icr_disabled)),1] <- "yellow"
colcol[c((length(icr_disabled)+1):(2*length(icr_disabled))),1] <- "green"

#Figure 3B
#=========================================================================================================
pdf("../Results/Revised_Figures/Pdfs/Common_TopMRs_Activity_ICR_Disabled_Figure_Revised(v2)_3B.pdf",height = 10, width=12, pointsize = 12)
heatmap.3(common_mrs_disabled_activity_matrix, Rowv = TRUE, Colv=TRUE, col = bluered(100), scale="none", main= "Median Activity of MRs common across all ICR Disabled Cancers",
          labRow = rownames(common_mrs_disabled_activity_matrix), labCol = colnames(common_mrs_disabled_activity_matrix), dendrogram = "col", 
          trace= "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.25, cexCol = 2, cellnote = ifelse(common_mrs_disabled_activity_matrix>0,"+","-"), notecex = 2, notecol = "black")
dev.off()

#Supplementary information about top common MRs for ICR Disabled cancer based on fold change values
A <- common_mrs_disabled_activity_matrix[common_mrs_disabled_list,c(1:length(icr_disabled))];
B <- common_mrs_disabled_activity_matrix[common_mrs_disabled_list,c((length(icr_disabled)+1):(2*length(icr_disabled)))]
wilcox_test_info <- perform_wilcox_test(A,B)
wilcox_test_info <- wilcox_test_info[order(wilcox_test_info$FC_Mean,decreasing=TRUE),]
wilcox_test_info <- revise_wilcox_test_info(wilcox_test_info)
write.table(wilcox_test_info,"../Results/Revised_Text_Results/Supplementary_Table_S5_Common_MRs_ICR_Disabled_Revised(v2).csv",row.names=T,col.names=T,quote=F)
top_10_common_mrs_disabled <- rownames(wilcox_test_info)[1:10]

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
  
  for (topmr in top_10_common_mrs_disabled)
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

p4 <- ggplot(data = disabled_activity_df, aes(x=TopMR, y=Activity_Value)) + 
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ Cancer, nrow=1, ncol=4) + xlab("Top Common MRs") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) + ylim(c(-0.75,0.75)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top 10 MRs common to ICR Disabled Cancers ") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))

ggsave(file="../Results/Revised_Figures/Pdfs/Common_TopMRs_Activity_Figure_Revised(v2)_3D.pdf",plot = p4, device = pdf(), width = 12, height=12, units = "in", dpi = 300)
dev.off()

###########################################################################################################################
#Load the images 3A and 3B
im1 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Common_TopMRs_Activity_ICR_Enabled_Figure_Revised(v2)_3A.svg")
p1 <- pictureGrob(im1)

im2 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Common_TopMRs_Activity_ICR_Disabled_Figure_Revised(v2)_3B.svg")
p2 <- pictureGrob(im2)

g_final <- ggarrange(p1,p3,p2,p4,nrow=2,ncol=2,labels=c("A","C","B","D"), widths = c(1.5,2.5))
ggsave(filename = paste0("../Results/Revised_Figures/Figure_Revised(v2)_3.pdf"),device = pdf(), plot = g_final, dpi = 300, width = 15, height= 11, units = "in" )
dev.off()

###########################################################################################################################
all_common_mrs <- intersect(common_mrs_enabled_list,common_mrs_disabled_list)