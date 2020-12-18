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
library(ggpubr)
library(grImport2)
library(matrixStats)
warnings("off")

#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
source("heatmap.3.R")

registerDoMC(20)

setwd('/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Old_Files/ICR_All_Info/scripts/')
setwd('.')

source('gene-reverse-network.R')
source('get_functions.R')

icr_enabled <- c("BLCA","BRCA","HNSC","LIHC","SARC","SKCM","STAD","UCEC")
icr_disabled <- c("LGG","KIRC","PAAD","UVM")
all_icr <- c(icr_enabled,icr_disabled)

#Make a set of density plot for the activity matrices of the 12 cancers of interest
#####################################################################################################
outputpath <- paste0("../Results/",all_icr[1]);
sample_name <- paste0(all_icr[1],"_Full_")
filename <- all_icr[1]
cols <- brewer.pal(12, "Paired")
line_types <- c(1,2,3,4,5,6,1,2,3,4,5,6)

#Load the activity matrix for BLCA
pdf("../Results/Revised_Figures/Pdfs/Normality_of_Activities.pdf",pointsize = 15, height=8, width=12)
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Activity_matrix_FGSEA.Rdata"))
plot(density(amat),xlab="Activity Values",ylab="Density",col=cols[1],type="l",lty=line_types[1],lwd=2,
     ylim=c(0,1),main="Activities of Transcriptional Regulators follows normal distribution", cex=2.0)
for (k in 2:length(all_icr))
{
  outputpath <- paste0("../Results/",all_icr[k]);
  sample_name <- paste0(all_icr[k],"_Full_")
  filename <- all_icr[k]
  #Load the activity matrix
  load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Activity_matrix_FGSEA.Rdata"))
  lines(density(amat),xlab="Activity Values",ylab="Density",col=cols[k],type="l",lty=line_types[k],lwd=2, cex=2.0)
}
legend("topright", legend=all_icr,col=cols, lty=line_types, cex=0.8)
dev.off()


#Make the plot of NES of the common MRs which are differentially active
output_df <- NULL
for (k in 1:length(all_icr))
{
  #Perform analysis for common MRs obtained from RGBM + FGSEA, RGBM + GSVA, RGBM + Viper, ARACNE + Viper
  #======================================================================================
  outputpath <- paste0("../Results/",all_icr[k]);
  sample_name <- paste0(all_icr[k],"_Full_")
  filename <- all_icr[k]
  
  #Load mechanistic network
  #======================================================================================
  load('../Data/Others/me_net_full.Rdata')
  
  #Load the gene expression matrix
  out <- loading_data(filename,M)
  D <- as.matrix(log2(t(out[[1]])+1))
  
  #Load the common MRs
  load(paste0(outputpath,"/Combined/",filename,"_Common_TopMRs.Rdata"))
  
  #Load NES information about all MRs identified by RGBM + FGSEA and intersect with common MRs
  load("../Results/BLCA/Adjacency_Matrix/BLCA_Full_TopMR_Info_FGSEA_BC_NES_1.Rdata")
  common_mr_info <- topmr_info[topmr_info$pathway %in% common_top_mrs,]
  common_mr_info$Cancer <- rep(filename,nrow(common_mr_info))
  
  #Load ICR high and ICR low indices
  load(paste0("../Data/",filename,"/",filename,"_ICR_cluster_assignment_k2-6.Rdata"))
  high_low_output <- get_high_low_indices(table_cluster_assignment,D)
  table_cluster_assignment <- high_low_output[[1]]
  high_indices <- high_low_output[[2]]
  low_indices <- high_low_output[[3]]
  D <- high_low_output[[4]]
  
  #Load the activity matrix
  load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Activity_matrix_FGSEA.Rdata"))
  max_amat <- max(amat)
  min_amat <- min(amat)
  for (i in 1:nrow(amat))
  {
     for (j in 1:ncol(amat))
     {
       if (amat[i,j]>0) {amat[i,j] <- amat[i,j]/max_amat}
       else if (amat[i,j]<0) {amat[i,j] <- amat[i,j]/abs(min_amat)}
     }
  }
  
  #common_mr_info$avg_icr_high <- rowSums(amat[common_mr_info$pathway,high_indices]/length(high_indices))
  #common_mr_info$avg_icr_low <- rowSums(amat[common_mr_info$pathway,low_indices]/length(low_indices))
  common_mr_info$median_icr_high <- rowMedians(amat[common_mr_info$pathway,high_indices])
  common_mr_info$median_icr_low <- rowMedians(amat[common_mr_info$pathway,low_indices])
  
  output_df <- rbind(output_df,cbind(common_mr_info$pathway,common_mr_info$padj,common_mr_info$NES,common_mr_info$Cancer,common_mr_info$median_icr_high,common_mr_info$median_icr_low))
}

output_df <- as.data.frame(output_df)
colnames(output_df) <- c("MR","Padj","NES","Cancer","Median_ICR_High","Median_ICR_Low")
output_df$MR <- as.character(as.vector(output_df$MR))
output_df$Padj <- as.numeric(as.vector(output_df$Padj))
output_df$NES <- as.numeric(as.vector(output_df$NES))
output_df$Cancer <- as.character(as.vector(output_df$Cancer))
output_df$Median_ICR_High <- as.numeric(as.vector(output_df$Median_ICR_High))
output_df$Median_ICR_Low <- as.numeric(as.vector(output_df$Median_ICR_Low))
write.table(output_df,"../Results/Revised_Text_Results/All_TF_Activity_Information_v2.csv",row.names = F, col.names = T, quote=F)


#Make a volcano plot highlighting the most differential TFs
#======================================================================================================================
output_df <- read.table("../Results/Revised_Text_Results/All_TF_Activity_Information_v2.csv",header = TRUE)
colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#a9a9a9', '#008080', '#e6beff')
output_df$MR <- as.character(as.vector(output_df$MR))
output_df$Cancer <- as.character(as.vector(output_df$Cancer))

temp_df <- output_df[-log10(output_df$Padj)>=4.0,]

p2 <- ggplot(output_df, aes(x = NES, y = -log10(Padj), 
                            shape = Cancer
)) + 
  geom_rect(data=NULL,aes(xmin=-4,xmax=0,ymin=-Inf,ymax=Inf),
            fill="green")+
  geom_rect(data=NULL,aes(xmin=0,xmax=4,ymin=-Inf,ymax=Inf),
            fill="yellow")+
  geom_point()+
  geom_point(data=output_df[output_df$Cancer==all_icr[12],], aes(x=NES, y=-log10(Padj), color=all_icr[12]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[1],], aes(x=NES, y=-log10(Padj), color=all_icr[1]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[2],], aes(x=NES, y=-log10(Padj), color=all_icr[2]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[3],], aes(x=NES, y=-log10(Padj), color=all_icr[3]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[4],], aes(x=NES, y=-log10(Padj), color=all_icr[4]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[5],], aes(x=NES, y=-log10(Padj), color=all_icr[5]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[6],], aes(x=NES, y=-log10(Padj), color=all_icr[6]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[7],], aes(x=NES, y=-log10(Padj), color=all_icr[7]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[8],], aes(x=NES, y=-log10(Padj), color=all_icr[8]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[9],], aes(x=NES, y=-log10(Padj), color=all_icr[9]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[10],], aes(x=NES, y=-log10(Padj), color=all_icr[10]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[11],], aes(x=NES, y=-log10(Padj), color=all_icr[11]))+
  xlab("Normalized Enrichment Scores (NES)") + ylab("-log10(Padj)") +
  scale_shape_manual(name = "Cancer", values = 1:nlevels(as.factor(output_df$Cancer))) +
  scale_color_manual(name = "Cancer", values =  colors) +
  geom_hline(yintercept=c(0,-log10(0.05),4),color=c("black","blue","red")) + geom_vline(xintercept = 0) + xlim(c(-4,4)) + ylim(c(0,6)) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background=element_rect(fill = "black"),
  #      panel.background = element_rect(fill = 'black'),legend.background = element_rect(fill = "black", color = NA),
  #      legend.key = element_rect(color = "gray", fill = "black"), legend.title = element_text(color = "white"), legend.text = element_text(color = "white"))+
  geom_text(aes(x=NES,y=-log10(Padj),label=MR),
            data=temp_df, 
            hjust = 0, vjust = 0, nudge_x = 0.01, angle = 45,
            size = 4, check_overlap = T) +
  ggtitle("Significance vs NES for Common MRs") + theme(text = element_text(size=16,color="black")) + theme(axis.text = element_text(size=14, color="black")) + theme(plot.title = element_text(hjust = 0.5,color="black")) + coord_fixed(ratio=0.7)
ggsave(file="../Results/Revised_Figures/Pdfs/Volcano_Plot_TF_Activity_Revised(v2)_2B.pdf",plot = p2, device = pdf(), width = 12, height=10, units = "in", dpi = 300)
dev.off()

supplementary_df <- output_df[,c(3,5,6)]
t_supp_df <- melt(supplementary_df,id.vars = "NES")
colnames(t_supp_df) <- c("NES","ICR","value")
t_supp_df$"NES_Eff" <- t_supp_df$NES>0
t_supp_df$ICR <- as.character(as.vector(t_supp_df$ICR))
t_supp_df[t_supp_df$ICR=="Median_ICR_High",]$ICR <- "ICR High"
t_supp_df[t_supp_df$ICR=="Median_ICR_Low",]$ICR <- "ICR Low"
t_supp_df[t_supp_df$NES>0,]$"NES_Eff" <- "NES>0"
t_supp_df[t_supp_df$NES<=0,]$"NES_Eff" <- "NES<0"

p_supp <- ggplot(t_supp_df, aes(x=NES_Eff,y=value, fill=ICR )) + geom_boxplot(position = position_dodge(0.9), outlier.alpha = 0.1 ) +
  scale_fill_manual(name="ICR Phenotype" ,values = c("yellow", "green")) + xlab("Categorized NES") + ylab("Median Common MR Activities")+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(c(-0.75,0.5)) +
  ggtitle("NES based classification of Common MRs") + theme(text = element_text(size=16)) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(x=1,y=0.5,label="****",col="red") + geom_text(x=2,y=0.5,label="****",col="red")
ggsave(file="../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_3a.jpg",plot=p_supp,device=jpeg(),width=12,height=10,units="in",dpi=300)
dev.off()

#Plot average activity of Common MRs coloring the MR per cancer
#============================================================================================================================
p3 <- ggplot(output_df, aes(x = Median_ICR_High, y = Median_ICR_Low, 
                            shape = Cancer,
                            size = -log10(Padj)
)) + 
  geom_point()+
  geom_point(data=output_df[output_df$Cancer==all_icr[12],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[12]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[1],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[1]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[2],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[2]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[3],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[3]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[4],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[4]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[5],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[5]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[6],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[6]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[7],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[7]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[8],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[8]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[9],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[9]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[10],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[10]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[11],], aes(x=Median_ICR_High, y=Median_ICR_Low, color=all_icr[11]))+
  xlab("Median Activity in ICR High samples") + ylab("Median Activity in ICR Low samples") +
  scale_shape_manual(name = "Cancer", values = 1:nlevels(as.factor(output_df$Cancer))) +
  scale_color_manual(name = "Cancer", values =  colors) +
  scale_size_continuous(range=c(0.1,3)) +
  geom_hline(yintercept=0,col="black") + geom_vline(xintercept = 0,col="black") + xlim(c(-0.4,0.4)) + ylim(c(-0.4,0.4)) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background=element_rect(fill = "black"),
  #      panel.background = element_rect(fill = 'black'),legend.background = element_rect(fill = "black", color = NA),
  #      legend.key = element_rect(color = "gray", fill = "black"), legend.title = element_text(color = "white"), legend.text = element_text(color = "white"))+
  geom_text(aes(x=Median_ICR_High,y=Median_ICR_Low,label=MR),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.01, angle = 45,
            size = 4, check_overlap = T,color="black") +
  ggtitle("Median Common MR Activity for ICR High vs ICR Low Phenotype") + theme(text = element_text(size=16,color="black")) + theme(axis.text = element_text(size=14, color="black")) + theme(plot.title = element_text(hjust = 0.5,color="black")) + coord_fixed(ratio=0.7)
ggsave(file="../Results/Revised_Figures/Pdfs/Avg_TF_Activity_Figure_Revised(v2)_2C.pdf",plot = p3, device = pdf(), width = 12, height=10, units = "in", dpi = 300)
dev.off()

#Make Revised Figure 2
#================================================================================================================
figure <- ggarrange(p2, p3,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
ggsave(filename = paste0("../Results/Revised_Figures/Figure_Revised(v2)_2.pdf"),device = pdf(), plot = figure, dpi = 300, width = 18, height= 8, units = "in" )
dev.off()
