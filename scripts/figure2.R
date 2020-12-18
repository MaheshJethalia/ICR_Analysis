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
library(grImport)
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
all_icr <- c(icr_enabled,icr_disabled)
output_df <- NULL

for (k in 1:length(all_icr))
{
  #Load the activity matrix of all TFs which satisfy minimum regulon criterion
  #======================================================================================
  outputpath <- paste0("Results/",all_icr[k]);
  sample_name <- paste0(all_icr[k],"_Full_")
  filename <- all_icr[k]
  
  #Load the mechanistic network
  load('Data/Others/me_net_full.Rdata')
  
  #Load the RNA-Seq matrix
  out <- loading_data(filename,M)
  D <- as.matrix(log2(t(out[[1]])+1))
  
  #Load the TFs with regulon 
  load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"regulon_FGSEA.Rdata"))
  
  #Load the ICR information
  load(paste0("Data/",filename,"/",filename,"_ICR_cluster_assignment_k2-6.Rdata"))
  tcga_sample_ids <- rownames(table_cluster_assignment)
  unmatched_tcga_sample_ids <- colnames(D);
  filter_samples <- which(unmatched_tcga_sample_ids %in% tcga_sample_ids);
  D <- D[,filter_samples];
  order_indices <- NULL
  for (i in 1:length(unmatched_tcga_sample_ids))
  {
    order_indices <- c(order_indices,which(tcga_sample_ids==unmatched_tcga_sample_ids[i]));
  }
  table_cluster_assignment <- table_cluster_assignment[order_indices,];
  high_indices <- which(table_cluster_assignment$HML_cluster=="ICR High")
  low_indices <- which(table_cluster_assignment$HML_cluster=="ICR Low")
  
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
  
  #Get difference in expression of target genes using hot vs cold phenotype
  pheno_t <- c()
  pv_t <- c()
  for (j in 1:nrow(D)){
    cat(j,"\n")
    high_sens <- mean(D[j,high_indices]);
    low_sens <- mean(D[j,low_indices]);
    pheno_t <- c(pheno_t,high_sens-low_sens);
    pv_t <- c(pv_t,wilcox.test(D[j,high_indices],D[j,low_indices],exact=F)$p.value)
  }
  padj_t = p.adjust(pv_t,method="fdr")
  names(pheno_t) <- rownames(D);
  sorted_pheno_t <- sort(pheno_t,decreasing = TRUE)
  
  #Perform gene set enrichment
  set.seed(123)
  fgseaRes <- fgsea(regulon_sets,sorted_pheno_t,minSize=10,maxSize=max(unlist(lapply(regulon_sets,length))), nperm=1000000)
  
  mra=as.data.frame(fgseaRes[order(pval)])
  mra$leadingEdge=as.character(mra$leadingEdge)
  if (ncol(D)>200)
  {
    topmr = mra$pathway[mra$padj<0.05 & abs(mra$NES)>1.75]
    topmr_info = mra[mra$padj<0.05 & abs(mra$NES)>1.75,]
  }
  else{
    topmr = mra$pathway[mra$padj<0.05 & abs(mra$NES)>1.0]
    topmr_info = mra[mra$padj<0.05 & abs(mra$NES)>1.0,]
  }
  save(topmr_info,file=paste0(outputpath,"/Adjacency_Matrix/",sample_name,"TopMR_Info_FGSEA_BC.Rdata"))
  mra$topmr_col <- rep(0,nrow(mra))
  mra[mra$pathway %in% topmr_info$pathway, ]$topmr_col <- 1
  mra$Cancer <- rep(all_icr[k],nrow(mra))
  mra$avg_icr_high <- rowSums(amat[mra$pathway,high_indices])/length(high_indices)
  mra$avg_icr_low <- rowSums(amat[mra$pathway,low_indices])/length(low_indices)
  
  output_df <- rbind(output_df,cbind(mra$pathway,mra$padj,mra$NES,mra$topmr_col,mra$Cancer,mra$avg_icr_high,mra$avg_icr_low))
}

output_df <- as.data.frame(output_df)
colnames(output_df) <- c("TF","Padj","NES","TopMR","Cancer","Avg_ICR_High","Avg_ICR_Low")
output_df$TF <- as.character(as.vector(output_df$TF))
output_df$Padj <- as.numeric(as.vector(output_df$Padj))
output_df$NES <- as.numeric(as.vector(output_df$NES))
output_df$TopMR <- as.numeric(as.vector(output_df$TopMR))
output_df$Cancer <- as.character(as.vector(output_df$Cancer))
output_df$Avg_ICR_High <- as.numeric(as.vector(output_df$Avg_ICR_High))
output_df$Avg_ICR_Low <- as.numeric(as.vector(output_df$Avg_ICR_Low))
write.table(output_df,"Results/Text_Results/All_TF_Activity_Information.csv",row.names = F, col.names = T, quote=F)


#Make heatmap for correlation based on common TF activites for 12 cancer for both ICR High and ICR Low cases
#===========================================================================================================
output_df <- read.table("Results/Text_Results/All_TF_Activity_Information.csv",header = TRUE)
colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#a9a9a9', '#008080', '#e6beff')
output_df$TF <- as.character(as.vector(output_df$TF))
output_df$Cancer <- as.character(as.vector(output_df$Cancer))
common_TFs <- output_df[output_df$Cancer==all_icr[1],]$TF
for (i in 2:length(all_icr))
{
  cancer <- all_icr[i]
  common_TFs <- intersect(common_TFs,output_df[output_df$Cancer==cancer,]$TF)
}
cor_matrix <- matrix(0,nrow=length(all_icr),ncol=length(all_icr))
cor_pval_matrix <- matrix(0,nrow=length(all_icr),ncol=length(all_icr))

for (k in 1:length(all_icr))
{
  for (l in 1:length(all_icr))
  {
    cancer1 <- all_icr[k]
    cancer2 <- all_icr[l]
    t1_df <- output_df[output_df$Cancer==cancer1,]
    t2_df <- output_df[output_df$Cancer==cancer2,]
    load(paste0("Results/",cancer1,"/Adjacency_Matrix/",cancer1,"_Full_Activity_matrix_FGSEA.Rdata"))
    amat[amat>0] <- amat[amat>0]/max(amat)
    amat[amat<0] <- amat[amat<0]/abs(min(amat))
    tf_cancer1_avg_activity <- rowMeans(amat[common_TFs,])
    
    load(paste0("Results/",cancer2,"/Adjacency_Matrix/",cancer2,"_Full_Activity_matrix_FGSEA.Rdata"))
    amat[amat>0] <- amat[amat>0]/max(amat)
    amat[amat<0] <- amat[amat<0]/abs(min(amat))
    tf_cancer2_avg_activity <- rowMeans(amat[common_TFs,])
    
    cor_matrix[k,l] <- cor(tf_cancer1_avg_activity,tf_cancer2_avg_activity,method = "pearson")
    test_info <- cor.test(tf_cancer1_avg_activity,tf_cancer2_avg_activity)
    cor_pval_matrix[k,l] <- test_info$p.value
    
  }
}
cor_pval_matrix <- matrix(p.adjust(as.vector(cor_pval_matrix),method="fdr"),nrow=length(all_icr),ncol=length(all_icr))
rownames(cor_matrix) <- all_icr
colnames(cor_matrix) <- all_icr
rownames(cor_pval_matrix) <- all_icr
colnames(cor_pval_matrix) <- all_icr
colcol <- matrix(0,nrow=length(all_icr),ncol=1)
colcol[c(1:length(icr_enabled)),1] <- "#FDB100"
colcol[c((length(icr_enabled)+1):length(all_icr)),1] <- "#660066"

pdf(paste0("Results/Figures/Correlation_Matrix_Activities_Figure2A.pdf"),width=14,height=14, pointsize = 16)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.5)
heatmap.3(cor_matrix, Rowv = TRUE, Colv="Rowv", col = bluered(100), scale="none", main= "Primary tumor similarity",
          labRow = rownames(cor_matrix), labCol = colnames(cor_matrix), dendrogram = "row", 
          key = TRUE, keysize=1.5, density.info = "none", KeyValueName = "Pearson Correlation", ColSideColors = colcol, ColSideColorsSize = 2,
          cellnote = ifelse(cor_pval_matrix<0.05, "*", NA), notecex = 2, notecol = "black", margins = c(6.5,6.5), useRaster = TRUE,
          cexRow = 2.5, cexCol = 2.5)
dev.off()

im1 <- readPicture("Results/Figures/Correlation_Matrix_Activities_Figure2A.svg")
p1 <- pictureGrob(im1)

#Make a volcano plot highlighting the most differential TFs
#======================================================================================================================
temp_df <- output_df[-log10(output_df$Padj)>=4.0,]

p2 <- ggplot(output_df, aes(x = NES, y = -log10(Padj), 
                            shape = Cancer
)) + 
  geom_rect(data=NULL,aes(xmin=-4,xmax=0,ymin=-Inf,ymax=Inf),
            fill="green")+
  geom_rect(data=NULL,aes(xmin=0,xmax=4,ymin=-Inf,ymax=Inf),
            fill="yellow")+
  geom_point()+
  geom_point(data=output_df[output_df$Cancer==all_icr[12] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[12]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[1] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[1]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[2] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[2]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[3] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[3]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[4] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[4]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[5] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[5]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[6] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[6]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[7] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[7]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[8] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[8]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[9] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[9]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[10] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[10]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[11] & output_df$TopMR==1,], aes(x=NES, y=-log10(Padj), color=all_icr[11]))+
  xlab("Normalized Enrichment Scores (NES)") + ylab("-log10(Padj)") +
  scale_shape_manual(name = "Cancer", values = 1:nlevels(as.factor(output_df$Cancer))) +
  scale_color_manual(name = "Cancer", values =  colors) +
  geom_hline(yintercept=c(0,2,4),color=c("black","blue","red")) + geom_vline(xintercept = 0) + xlim(c(-4,4)) + ylim(c(0,6)) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background=element_rect(fill = "black"),
  #      panel.background = element_rect(fill = 'black'),legend.background = element_rect(fill = "black", color = NA),
  #      legend.key = element_rect(color = "gray", fill = "black"), legend.title = element_text(color = "white"), legend.text = element_text(color = "white"))+
  geom_text(aes(x=NES,y=-log10(Padj),label=TF),
            data=temp_df, 
            hjust = 0, vjust = 0, nudge_x = 0.01, angle = 45,
            size = 4, check_overlap = T) +
  ggtitle("Significance vs NES for TF activities") + theme(text = element_text(size=16,color="black")) + theme(axis.text = element_text(size=14, color="black")) + theme(plot.title = element_text(hjust = 0.5,color="black")) + coord_fixed(ratio=0.7)
ggsave(file="Results/Figures/Volcano_Plot_TF_Activity_2B.pdf",plot = p2, device = pdf(), width = 12, height=10, units = "in", dpi = 300)
dev.off()


supplementary_df <- output_df[,c(3,6,7)]
t_supp_df <- melt(supplementary_df,id.vars = "NES")
colnames(t_supp_df) <- c("NES","ICR","value")
t_supp_df$"NES_Eff" <- t_supp_df$NES>0
t_supp_df$ICR <- as.character(as.vector(t_supp_df$ICR))
t_supp_df[t_supp_df$ICR=="Avg_ICR_High",]$ICR <- "ICR High"
t_supp_df[t_supp_df$ICR=="Avg_ICR_Low",]$ICR <- "ICR Low"
t_supp_df[t_supp_df$NES>0,]$"NES_Eff" <- "NES>0"
t_supp_df[t_supp_df$NES<=0,]$"NES_Eff" <- "NES<0"

p_supp <- ggplot(t_supp_df, aes(x=NES_Eff,y=value, fill=ICR )) + geom_boxplot(position = position_dodge(0.9), outlier.alpha = 0.1 ) +
          scale_fill_manual(name="ICR Phenotype" ,values = c("yellow", "green")) + xlab("Categorized NES") + ylab("Avg TF Activities")+
          theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(c(-0.75,0.5)) +
          ggtitle("NES based classification of TFs") + theme(text = element_text(size=16)) + theme(plot.title = element_text(hjust = 0.5)) +
          geom_text(x=1,y=0.5,label="****",col="red") + geom_text(x=2,y=0.5,label="****",col="red")
ggsave(file="Results/Figures/Supp_Figure2.jpg",plot=p_supp,device=jpeg(),width=12,height=10,units="in",dpi=300)
dev.off()
          
#Plot average activity of TFs coloring the MR per cancer
#============================================================================================================================
p3 <- ggplot(output_df, aes(x = Avg_ICR_High, y = Avg_ICR_Low, 
                            shape = Cancer,
                            size = -log10(Padj)
)) + 
  geom_point()+
  geom_point(data=output_df[output_df$Cancer==all_icr[12] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[12]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[1] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[1]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[2] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[2]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[3] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[3]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[4] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[4]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[5] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[5]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[6] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[6]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[7] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[7]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[8] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[8]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[9] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[9]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[10] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[10]))+
  geom_point(data=output_df[output_df$Cancer==all_icr[11] & output_df$TopMR==1,], aes(x=Avg_ICR_High, y=Avg_ICR_Low, color=all_icr[11]))+
  xlab("Average Activity in ICR High samples") + ylab("Average Activity in ICR Low samples") +
  scale_shape_manual(name = "Cancer", values = 1:nlevels(as.factor(output_df$Cancer))) +
  scale_color_manual(name = "Cancer", values =  colors) +
  scale_size_continuous(range=c(0.1,3)) +
  geom_hline(yintercept=0,col="black") + geom_vline(xintercept = 0,col="black") + xlim(c(-0.5,0.5)) + ylim(c(-0.3,0.3)) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background=element_rect(fill = "black"),
  #      panel.background = element_rect(fill = 'black'),legend.background = element_rect(fill = "black", color = NA),
  #      legend.key = element_rect(color = "gray", fill = "black"), legend.title = element_text(color = "white"), legend.text = element_text(color = "white"))+
  geom_text(aes(x=Avg_ICR_High,y=Avg_ICR_Low,label=TF),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.01, angle = 45,
            size = 4, check_overlap = T,color="black") +
  ggtitle("Average TF Activity for ICR High vs ICR Low Phenotype") + theme(text = element_text(size=16,color="black")) + theme(axis.text = element_text(size=14, color="black")) + theme(plot.title = element_text(hjust = 0.5,color="black")) + coord_fixed(ratio=0.7)
ggsave(file="Results/Figures/Pdfs/Avg_TF_Activity_Figure2C.pdf",plot = p3, device = pdf(), width = 12, height=10, units = "in", dpi = 300)
dev.off()


#Get TopMR and subnetwork information for BRCA and visuallize the same
#=================================================================================================================
filename <- "BRCA"
outputpath <- paste0("Results/",filename);
sample_name <- paste0(filename,"_Full_")
net <- read.table(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv"),header=TRUE,sep=" ")

load("Data/Others/TF_Target_Info.Rdata")
tfs <- tf_gene_names;
targets <- target_gene_names;
rownames(net) <- tfs;
colnames(net) <- targets;

load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"TopMR_Info_FGSEA_BC.Rdata"))
topmrs <- topmr_info$pathway

topmr_subnet <- net[topmrs,]
topmr_edgelist <- NULL
for (i in 1:nrow(topmr_subnet))
{
  topmr <- topmrs[i]
  targets <- colnames(topmr_subnet[topmr,which(topmr_subnet[topmr,]>0)])
  temp <- cbind(rep(topmr,length(targets)),targets,rep(1,length(targets)))
  topmr_edgelist <- rbind(topmr_edgelist,temp)
}

topmr_edgelist <- as.data.frame(topmr_edgelist)
colnames(topmr_edgelist) <- c("Source","Target","Weight")
topmr_edgelist$Source <- as.character(as.vector(topmr_edgelist$Source))
topmr_edgelist$Target <- as.character(as.vector(topmr_edgelist$Target))
topmr_edgelist$Weight <- as.numeric(as.vector(topmr_edgelist$Weight))
write.table(topmr_edgelist,"Results/BRCA/Adjacency_Matrix/BRCA_TopMR_Subnetwork.csv",row.names=F,col.names=T,quote=F)

#Make the node table
all_targets <- colnames(net)
all_targets <- union(topmrs,setdiff(all_targets,topmrs))
node_info <- cbind(all_targets,as.numeric(all_targets %in% topmrs),c(as.numeric(output_df[output_df$Cancer==filename & output_df$TopMR==1, ]$NES>0)+1,
                                                                     rep(0,length(all_targets)-length(topmrs))))
node_info <- as.data.frame(node_info)
colnames(node_info) <- c("Node","TopMR","ICR")
node_info$Node <- as.character(as.vector(node_info$Node))
node_info$TopMR <- as.numeric(as.vector(node_info$TopMR))
node_info$ICR <- as.numeric(as.vector(node_info$ICR))
write.table(node_info,"Results/BRCA/Adjacency_Matrix/BRCA_TopMR_Nodeinfo.csv",row.names=F,col.names=T,quote=F)

im4 <- readPicture(paste0("Results/Figures/BRCA_TopMR_subnetwork_Figure_2D.svg"))
p4 <- pictureGrob(im4)

figure <- ggarrange(p1, p2, p3, p4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
ggsave(filename = paste0("Results/Figures/Figure2.pdf"),device = pdf(), plot = figure, dpi = 300, width = 18, height= 12, units = "in" )
dev.off()
