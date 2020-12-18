library(data.table)
library(ggplot2)
library(doMC)
library(heatmap3)
library(gplots)
library(igraph)
library(yaGST)
library(parallel)
#BiocManager::install("fgsea")
require(fgsea)
require(limma)
registerDoMC(20)

setwd('/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Old_Files/ICR_All_Info/scripts/')
setwd('.')

source('gene-reverse-network.R')
source('get_functions.R')

#Load mechanistic network
#======================================================================================
load('../Data/Others/me_net_full.Rdata')

#Load the RNASeq data
#=======================================================================================
filename <- "KIRC"
outputpath <- paste0("../Results/",filename)
sample_name <- paste0(filename,"_Full_");
out <- loading_data(filename,M)
D <- as.matrix(log2(t(out[[1]])+1))

#Get regulatory network
#=======================================================================================
net <- read.table(gzfile(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv.gz")),header=TRUE,sep=" ")

#Get the gene names and associate it with appropriate variables
#====================================================================================
load("../Data/Others/TF_Target_Info.Rdata")
tfs <- tf_gene_names;
targets <- target_gene_names;
rownames(net) <- tfs;
colnames(net) <- targets;

#Make average regulon size to less than 100
#====================================================================================
while((sum(net>0)/length(tfs))>100)
{
  median_weight <- median(net[net>0]);
  net[net<median_weight] <- 0;
}


#Get phenotype only for samples which have phenotype information
#=====================================================================================
load(paste0("../Data/",filename,"/",filename,"_ICR_cluster_assignment_k2-6.Rdata"))
high_low_output <- get_high_low_indices(table_cluster_assignment,D)
table_cluster_assignment <- high_low_output[[1]]
high_indices <- high_low_output[[2]]
low_indices <- high_low_output[[3]]
D <- high_low_output[[4]]

#Get the tumor purity information for each tumor sample
#######================================================================================
tumor_purity_info <- fread("../Data/Others/NIHMS958047-supplement-1.csv",header=TRUE,sep="\t")
tumor_purity_info <- as.data.frame(tumor_purity_info)
tumor_purity_info$Type <- as.character(as.vector(tumor_purity_info$Type))
tumor_purity_info$Sample <- as.character(as.vector(tumor_purity_info$Sample))

#Get the subset of tumor purity information for a particular cancer
tumor_purity_subset <- tumor_purity_info[tumor_purity_info$Type==filename,]
tumor_purity_vec <- rep(1.0,dim(D)[2])

for (i in 1:nrow(tumor_purity_subset))
{
  sample <- tumor_purity_subset$Sample[i]
  purity_level <- tumor_purity_subset$Purity[i]
  if (!is.na(purity_level))
  {
    ids <- which(grepl(sample,colnames(D),fixed=TRUE))
    tumor_purity_vec[ids] <- purity_level
  }
}
tumor_purity_matrix <- matrix(tumor_purity_vec, nrow=dim(D)[1], ncol=length(tumor_purity_vec), byrow=TRUE)
colnames(tumor_purity_matrix) <- colnames(D)
rownames(tumor_purity_matrix) <- rownames(D)

#Create activity matrix
#=====================================================================================
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Activity_matrix_FGSEA.Rdata"))

#Get the ordered activity matrix for significant MRs
#==============================================================================================
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"TopMR_Info_FGSEA_BC_NES_1.Rdata"))
topmr <- as.character(as.vector(topmr_info$pathway))
amat.mr <- amat[topmr,];

#Perform limma without tumor purity as covariate
amat_signal <- amat.mr[,c(high_indices,low_indices)];
groups <- c(rep("ICR_High",length(high_indices)),rep("ICR_Low",length(low_indices)))
fit_without_tp <- eBayes(lmFit(amat_signal, model.matrix(~groups)))
limma_without_tp <- topTable(fit_without_tp, coef=2, number = Inf)

#Perform limma taking into account the effect of tumor purity
tumor_purity_signal <- tumor_purity_vec[c(high_indices,low_indices)]
fit_with_tp <- eBayes(lmFit(amat_signal, model.matrix(~tumor_purity_signal+groups)))
limma_with_tp <- topTable(fit_with_tp, coef=3, number=Inf)
#write.table(limma_with_tp,file=paste0(outputpath,"/Adjacency_Matrix/",sample_name,"WilCox_Test_FGSEA_With_Tumor_Purity_MR_Result_NES_1.csv"),row.names=T,col.names=T);

#Compare the results for the case without and with tumor purity
temp2 <- limma_without_tp[order(rownames(limma_without_tp)),]
temp1 <- limma_with_tp[order(rownames(limma_with_tp)),]
all_diff_tfs <- union(rownames(temp2[temp2$adj.P.Val<0.05,]),rownames(temp1[temp1$adj.P.Val<0.05,]))
temp3 <- temp2[all_diff_tfs,]
temp4 <- temp1[all_diff_tfs,]
x <- -log10(temp3$adj.P.Val)
y <- -log10(temp4$adj.P.Val)
new_df <- data.frame(x=x,y=y)
reg<-lm(y~x,data=new_df)

pdf(paste0(outputpath,"/Images/Comparison_of_MRs_for_",filename,"_tumor_purity_NES_1.pdf"),height=7,width=10)
plot(new_df$x,new_df$y,xlab="-log10(P.adj) of Top MRs without tumor purity",
     ylab="-log10(P.adj) of Top MRs with tumor purity", 
     main=paste0("Comparison of differential activity for ",filename," without or with tumor purity"),
     cex.main=1.5, cex.lab=1.5, cex.axis=1.0)
abline(h=-log10(0.05),col="red",lwd=2,lty=2)
abline(v=-log10(0.05),col="red",lwd=2,lty=2)
abline(reg, col="blue",lwd=2,lty=1)
dev.off()

library(VennDiagram)
significant_mrs_with_tumor_purity <- rownames(temp1[temp1$adj.P.Val<=0.05,])
significant_mrs_without_tumor_purity <- rownames(temp2[temp2$adj.P.Val<=0.05,])
venn.diagram(x = list(significant_mrs_with_tumor_purity,significant_mrs_without_tumor_purity),
             category.names = c("With TP", "Without TP"),
             height = 1200,
             width = 1200,
             resolution = 300,
             imagetype = "tiff",
             fill = c("#9297FF", "#FFD38F"),
             fontfamily = "Arial",
             alpha = c(0.5, 0.5),
             main = "Comparison of Top MRs",
             filename = paste0(outputpath,"/Images/Intersection_of_MRs_for_",filename,"_tumor_purity_NES_1.tiff"),
             output=TRUE
)