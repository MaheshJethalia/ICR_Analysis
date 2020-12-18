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
library(gridExtra)
library(grImport2)
library(ggpubr)

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

registerDoMC(20)

setwd('/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Old_Files/ICR_All_Info/scripts/')
setwd('.')

source('gene-reverse-network.R')
source('get_functions.R')

#Load mechanistic network
#======================================================================================
load('../Data/Others/me_net_full.Rdata')

#Load the RNASeq data for BLCA
#=======================================================================================
filename <- "BRCA"
outputpath <- paste0("../Results/",filename)
sample_name <- paste0(filename,"_Full_");
out <- loading_data(filename,M)
D <- as.matrix(log2(t(out[[1]])+1))

#Filter to get the samples with the ICR phenotype information
#======================================================================================
load(paste0("../Data/",filename,"/",filename,"_ICR_cluster_assignment_k2-6.Rdata"))
tcga_sample_ids <- rownames(table_cluster_assignment)
unmatched_tcga_sample_ids <- colnames(D);
filter_samples <- which(unmatched_tcga_sample_ids %in% tcga_sample_ids);
D <- D[,filter_samples];
jpeg(filename="../Results/Revised_Figures/Svgs_Jpgs/Supp_Figure_Revised(v2)_2.jpg",width=900, height=480, units="px", pointsize=12)
boxplot(D[,1:100],ylab="Quantile-Normalized Counts",main="Quantile-Normalized Primary Tumor Samples")
dev.off()

order_indices <- NULL
for (i in 1:length(unmatched_tcga_sample_ids))
{
  order_indices <- c(order_indices,which(tcga_sample_ids==unmatched_tcga_sample_ids[i]));
}
table_cluster_assignment <- table_cluster_assignment[order_indices,];
high_indices <- which(table_cluster_assignment$HML_cluster=="ICR High")
low_indices <- which(table_cluster_assignment$HML_cluster=="ICR Low")
D1 <- D[,c(high_indices,low_indices)]
x <- apply(D1,1,IQR)
D2 <- D1[x>quantile(x,0.95),]

##Perform hierarchical clustering to get clusters on both rows and columns
hc.cols <- hclust(dist(t(D2)),method="ward.D2")
hc.rows <- hclust(dist(D2),method="ward.D2")

pdf(paste0("../Results/Revised_Figures/Pdfs/RNA_Seq_Figure_Revised(v2)_1A.pdf"),width=12,height=10,pointsize = 16)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="white", cex.main=2.0)
heatmap.3(D2[hc.rows$order,hc.cols$order], Rowv = NULL, Colv = NULL, col = bluered(50), scale="none", main= "RNA-Seq Expression Data",
          xlab = "TCGA Samples", ylab = "Genes ", labRow = FALSE, labCol = FALSE, dendrogram = "none", cexRow = 6, cexCol = 6,
          key = TRUE, density.info = "none", KeyValueName = "log2(Normalized Expression)",
          margins = c(2,2), useRaster = TRUE)
dev.off()

#Load cluster RNA-Seq matrix image
im1 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/RNA_Seq_Figure_Revised(v2)_1A.svg")
g1 <- pictureGrob(im1)
###################################################################################################################

#Get regulatory network
#=======================================================================================
net <- read.table(gzfile(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv.gz")),header=TRUE,sep=" ")

#Get the gene names and associate it with appropriate variables
#====================================================================================
load("../Data/Others/TF_Target_Info.Rdata")
tfs <- tf_gene_names;
targets <- target_gene_names;
rownames(net) <- tfs;
colnames(net) <- targets

#Load inferred GRN for BLCA
im2 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/BLCA_GRN_Figure_Revised(v2)_1B.svg")
g2 <- pictureGrob(im2)

###################################################################################################################

#Load the TF-target signed regulon
#===================================================================================
im4 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/Toy_TF_Regulon_Figure_Revised(v2)_1C.svg")
g4 <- pictureGrob(im4)

#######################################################################################################################
#Load the activity matrix for all TFs with regulon size >= 10
#==================================================================================
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Activity_matrix_FGSEA.Rdata"))
amat.pheno <- amat[,c(high_indices,low_indices)]
hc.cols <- hclust(dist(t(amat.pheno)),method="ward.D2")
hc.rows <- hclust(dist(amat.pheno),method="ward.D2")

min_activity <- min(amat.pheno)
max_activity <- max(amat.pheno)
for (i in 1:ncol(amat.pheno))
{
  for (j in 1:nrow(amat.pheno))
  {
    if (amat.pheno[j,i]>0)
    {
      amat.pheno[j,i] <- amat.pheno[j,i]/max_activity
    } else{
      amat.pheno[j,i] <- amat.pheno[j,i]/abs(min_activity)
    }
  }
  
}
hc.cols <- hclust(dist(t(amat.pheno)),method="ward.D2")
hc.rows <- hclust(dist(amat.pheno),method="ward.D2")

pdf("../Results/Revised_Figures/Pdfs/BLCA_TF_Activity_Figure_Revised(v2)_1D.pdf",width=12,height=12,pointsize=16)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0, cex = 2.5)
heatmap.3(amat.pheno[hc.rows$order,hc.cols$order], Rowv = NULL, Colv = NULL, col = bluered(50), scale="none", main= "TF Regulon Activity Matrix",
          xlab = "TCGA Samples", ylab = "TFs", labRow = FALSE, labCol = FALSE, dendrogram = "none", cexRow=6, cexCol=6,
          key = TRUE, density.info = "none", KeyValueName = "Activity Value",
          margins = c(2,2), useRaster = TRUE)
dev.off()

#Load the activity matrix for TFs with regulon >=10
#===================================================================================
im5 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/BLCA_TF_Activity_Figure_Revised(v2)_1D.svg")
g5 <- pictureGrob(im5)

###########################################################################################################################

#Get difference in expression of target genes using hot vs cold phenotype
#===================================================================================
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
#========================================================================================
require(fgsea)
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"regulon_FGSEA.Rdata"))
fgseaRes <- fgsea(regulon_sets,sorted_pheno_t,minSize=10,maxSize=max(unlist(lapply(regulon_sets,length))), nperm=1000000)
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#Get plot of top MRs and select top MRs
#========================================================================================
pdf("../Results/Revised_Figures/Pdfs/BLCA_GeneRanks_FGSEA_Figure_Revised(v2)_1E.pdf",pointsize=18,height = 10, width = 10)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0, cex = 2.5)
plotGseaTable(regulon_sets[topPathways], sorted_pheno_t, fgseaRes, gseaParam = 0.5)
dev.off()

im6 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/BLCA_GeneRanks_FGSEA_Figure_Revised(v2)_1E.svg")
g6 <- pictureGrob(im6)

######################################################################################################################
#Get Top MR make their activity matrix
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"TopMR_Info_FGSEA_BC.Rdata"))
new_mr_activity_matrix <- as.matrix(amat[topmr_info$pathway,c(high_indices,low_indices)])
rownames(new_mr_activity_matrix)
colnames(new_mr_activity_matrix) <- NULL
min_activity <- min(new_mr_activity_matrix)
max_activity <- max(new_mr_activity_matrix)
for (i in 1:ncol(new_mr_activity_matrix))
{
  for (j in 1:nrow(new_mr_activity_matrix))
  {
    if (new_mr_activity_matrix[j,i]>0)
    {
      new_mr_activity_matrix[j,i] <- new_mr_activity_matrix[j,i]/max_activity
    } else{
      new_mr_activity_matrix[j,i] <- new_mr_activity_matrix[j,i]/abs(min_activity)
    }
  }
  
}

colcol <- matrix(0,nrow=ncol(new_mr_activity_matrix),ncol=1)
#rep("white",ncol(new_mr_activity_matrix));
high_ids <- c(1:length(high_indices));
low_ids <- c((length(high_indices)+1):length(colcol))
colcol[high_ids,1] <- "yellow"
colcol[low_ids,1] <- "green"
colnames(colcol) <- "ICR Phenotype"
hc_high.cols <- hclust(dist(t(new_mr_activity_matrix[,high_ids])),method='ward.D')
hc_low.cols <- hclust(dist(t(new_mr_activity_matrix[,low_ids])),method='ward.D')
hc.cols <- c(hc_high.cols$order,length(high_ids)+hc_low.cols$order);
hc.rows <- hclust(dist(new_mr_activity_matrix),method='ward.D')

#Make the plot of activity matrix clustered by ICR High vs ICR Low
#=================================================================================================
pdf("../Results/Revised_Figures/Pdfs/BLCA_TopMR_Activity_Figure_Revised(v2)_1F.pdf",height=10,width=12,pointsize=16)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0, cex = 2.5)
heatmap.3(new_mr_activity_matrix[hc.rows$order,hc.cols], Rowv = NULL, Colv = NULL, col = bluered(50), scale="none", main= "Master Regulator Activity Matrix",
          xlab = "TCGA Samples", ylab = "Top MRs", labRow = FALSE, labCol = FALSE, dendrogram = "none",
          key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize=2,cexCol=2,
          margins = c(2,2), useRaster = TRUE)
dev.off()

im7 <- readPicture("../Results/Revised_Figures/Svgs_Jpgs/BLCA_TopMR_Activity_Figure_Revised(v2)_1F.svg")
g7 <- pictureGrob(im7)

#g_final <- grid.arrange( grobs = list(g1,g2,g4,g5,g6,g7),  widths = c(1, 1.25, 1.5), layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))
g_final <- ggarrange(g1,g2,g4,g5,g6,g7, 
                     labels = c("A","B","C","D","E","F"),
                     nrow=2,ncol=3)


ggsave(filename = paste0("../Results/Revised_Figures/Figure_Revised(v2)_1.pdf"),device = pdf(), plot = g_final, dpi = 300, width = 12, height= 10.5, units = "in" )
dev.off()

