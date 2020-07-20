library(data.table)
library(ggplot2)
library(gplots)
library(Matrix)
library(erer)
library(doMC)
library(heatmap3)
library(mgcv)
library(devtools)

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

setwd('../ICR_All_Info/')

all_cold_mrs_df <- read.table("Results/Figures/Supp_Figure6.csv",header=TRUE,sep=",")
all_cold_mrs_df <- all_cold_mrs_df[all_cold_mrs_df$Avg_NES_All_Cancers<0,]
all_cold_mrs <- as.character(as.vector(all_cold_mrs_df[,1]))

all_hot_mrs_df <- read.table("Results/Figures/Supp_Figure5.csv",header=TRUE,sep=",")
all_hot_mrs_df <- all_hot_mrs_df[all_hot_mrs_df$Avg_NES_All_Cancers>0,]
all_hot_mrs <- as.character(as.vector(all_hot_mrs_df[,1]))
all_mrs <- c(all_hot_mrs,all_cold_mrs)

RGBM_path <- c("Data/");

icr_disabled_cancers <- c("LGG","PAAD","UVM","KIRC")
icr_enabled_cancers <- c("BLCA","BRCA","HNSC","LIHC","SARC","STAD","SKCM","UCEC")
icr_neutral_cancers <- setdiff(setdiff(setdiff(list.files(RGBM_path),icr_disabled_cancers),icr_enabled_cancers),c("Others"))

get_activity_matrix_info <- function(counter,sample_ids,cancers,activity_matrix,cancer_info_vector,icr_pheno_info_vector,icr_type_info_vector,icr_type)
{
  for (i in 1:length(cancers))
  {
    cancer_type <- cancers[i]
    load(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
    
    high_indices_table <- read.table(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_high_indices.csv"),header=TRUE)
    low_indices_table <- read.table(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_low_indices.csv"),header=TRUE)
    amat_to_use <- amat[,c(high_indices_table$x,low_indices_table$x)]
    amat_to_use[amat_to_use>0] <- amat_to_use[amat_to_use>0]/max(amat_to_use)
    amat_to_use[amat_to_use<0] <- amat_to_use[amat_to_use<0]/abs(min(amat_to_use))
    samples <- dim(amat_to_use)[2]
    high_index_last <- nrow(high_indices_table)
    low_index_last <- ncol(amat_to_use)
    
    cancer_info_vector <- c(cancer_info_vector,rep(cancer_type,samples))
    icr_pheno_info_vector <- c(icr_pheno_info_vector,rep("ICR High",nrow(high_indices_table)),
                               rep("ICR Low",nrow(low_indices_table)))
    icr_type_info_vector <- c(icr_type_info_vector,rep(icr_type,samples))
    sample_ids <- c(sample_ids,colnames(amat_to_use))
    for (j in 1:length(all_mrs))
    {
      mr <- all_mrs[j]
      if (mr %in% rownames(amat_to_use))
      {
        activity_matrix[mr,counter+c(1:high_index_last)] <- amat_to_use[mr,c(1:high_index_last)]
        activity_matrix[mr,counter+c((high_index_last+1):low_index_last)] <- amat_to_use[mr,c((high_index_last+1):low_index_last)]
      }
    }
    counter <- counter+samples
  }
  
  output_list <- list()
  output_list[[1]] <- counter
  output_list[[2]] <- cancer_info_vector
  output_list[[3]] <- icr_pheno_info_vector
  output_list[[4]] <- icr_type_info_vector
  output_list[[5]] <- activity_matrix
  output_list[[6]] <- sample_ids
  return(output_list)
}


cancers <- c(icr_disabled_cancers,icr_enabled_cancers,icr_neutral_cancers)
cancer_samples <- 0
for (i in 1:length(cancers))
{
  cancer_type <- cancers[i]
  high_indices_table <- read.table(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_low_indices.csv"),header=TRUE)
  cancer_samples <- cancer_samples+nrow(high_indices_table)+nrow(low_indices_table)
}

Activity_Matrix <- Matrix(0,nrow=length(all_mrs),ncol=cancer_samples)
rownames(Activity_Matrix) <- all_mrs
counter=0;
cancer_info_vector <- NULL
icr_pheno_info_vector <- NULL
icr_type_info_vector <- NULL
sample_ids <- NULL

output <- get_activity_matrix_info(counter,sample_ids,icr_disabled_cancers,Activity_Matrix,cancer_info_vector,icr_pheno_info_vector,icr_type_info_vector,"ICR Disabled")
counter <- output[[1]]
cancer_info_vector <- output[[2]]
icr_pheno_info_vector <- output[[3]]
icr_type_info_vector <- output[[4]]
Activity_Matrix <- output[[5]]
sample_ids <- output[[6]]

output <- get_activity_matrix_info(counter,sample_ids,icr_enabled_cancers,Activity_Matrix,cancer_info_vector,icr_pheno_info_vector,icr_type_info_vector,"ICR Enabled")
counter <- output[[1]]
cancer_info_vector <- output[[2]]
icr_pheno_info_vector <- output[[3]]
icr_type_info_vector <- output[[4]]
Activity_Matrix <- output[[5]]
sample_ids <- output[[6]]

output <- get_activity_matrix_info(counter,sample_ids,icr_neutral_cancers,Activity_Matrix,cancer_info_vector,icr_pheno_info_vector,icr_type_info_vector,"ICR Neutral")
counter <- output[[1]]
cancer_info_vector <- output[[2]]
icr_pheno_info_vector <- output[[3]]
icr_type_info_vector <- output[[4]]
Activity_Matrix <- output[[5]]
sample_ids <- output[[6]]
colnames(Activity_Matrix) <- sample_ids
save(Activity_Matrix,all_cold_mrs,all_hot_mrs,cancer_info_vector,icr_type_info_vector,file="Results/PanCancer/PanCancer_Activity_Info.Rdata")

#Make heatmap using heatmap.3 function
n <- length(unique(cancer_info_vector))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(43)
col=sample(color, n)
Cancers <- col[as.factor(cancer_info_vector)]
icr_colors <- c("yellow","green")
ICR_Clusters <- icr_colors[as.factor(icr_pheno_info_vector)]
icr_type_colors <- c("#660066","#FDB100","#D3D3D3")
ICR_EorD <- icr_type_colors[as.factor(icr_type_info_vector)]

clab=cbind(Cancers,ICR_Clusters,ICR_EorD)
rlab = t(c(rep("yellow",length(all_hot_mrs)),rep("green",length(all_cold_mrs))))
colnames(clab) <- c("Cancers","ICR Phenotype","ICR Cancer Clusters")
rownames(rlab) <- "MRs"

#Define custom dist and hclust functions for use with heatmaps
Activity_Matrix <- as.matrix(Activity_Matrix)
#2nd column is ICR High vs ICR Low
hc.cols1 <- order(clab[,2],decreasing = T)
clab1 <- clab[hc.cols1,]
Activity_Matrix <- Activity_Matrix[,hc.cols1]
newhc.cols <- NULL
#3rd column is ICR Enabled, ICR Disabled and ICR Neutral clusters
for (i in 1:length(unique(clab[,2])))
{
  color_id <- unique(clab[,2])[i]
  newhc.cols <- c(newhc.cols,length(newhc.cols)+order(clab1[clab1[,2]==color_id,3],decreasing = T))
}
Activity_Matrix <- Activity_Matrix[,newhc.cols]
clab2 <- clab1[newhc.cols,]
lower_val <- quantile(Activity_Matrix,probs=c(0.05))
upper_val <- quantile(Activity_Matrix,probs=c(0.95))
Activity_Matrix[Activity_Matrix<lower_val] <- lower_val
Activity_Matrix[Activity_Matrix>upper_val] <- upper_val

.merge_hclust <- function(hclist) {
  #-- Merge
  d <- as.dendrogram(hclist[[1]])
  for (i in 2:length(hclist)) {
    d <- merge(d, as.dendrogram(hclist[[i]]))
  }
  return(as.hclust(d))
}

all_icr_high_indices <- which(clab[,2]==icr_colors[1])
all_icr_low_indices <- which(clab[,2]==icr_colors[2])
hc.cols.icr_high <- hclust(dist(t(Activity_Matrix[,all_icr_high_indices])),method='ward.D2')
hc.cols.icr_low <- hclust(dist(t(Activity_Matrix[,all_icr_low_indices])),method='ward.D2')
hc.cols <- c(hc.cols.icr_high$order,(length(all_icr_high_indices)+hc.cols.icr_low$order))
hc_col_list <- list(hc.cols.icr_high,hc.cols.icr_low)
hc_col <- .merge_hclust(hc_col_list)
#clab3 <- clab2[hc.cols,]
#Activity_Matrix <- Activity_Matrix[,hc.cols]

all_hot_mr_indices <- which(rlab[1,]==icr_colors[1])
all_cold_mr_indices <- which(rlab[1,]==icr_colors[2])
hc.rows.icr_high <- hclust(dist(Activity_Matrix[all_hot_mr_indices,]),method='ward.D2')
hc.rows.icr_low <- hclust(dist(Activity_Matrix[all_cold_mr_indices,]),method='ward.D2')
hc.rows <- c(hc.rows.icr_high$order,(length(all_hot_mr_indices)+hc.rows.icr_low$order))
hc_row_list <- list(hc.rows.icr_high,hc.rows.icr_low)
hc_row <- .merge_hclust(hc_row_list)
#rlab2 <- t(rlab[,hc.rows])
#rownames(rlab2) <- "ICR MRs"
#Activity_Matrix <- Activity_Matrix[hc.rows,]


#Main Plotting Function
pdf("Results/Figures/Heatmap_All_MRs_PanCancer_Figure.pdf",height = 13, width=15, pointsize = 14)
par(bg="black")
par(fg="white",col.axis="white",col.main="white",col.lab="white", cex.main=2.0)
heatmap.3(Activity_Matrix, scale="none", dendrogram="both", margins=c(4,20),
          Rowv=as.dendrogram(hc_row), Colv=as.dendrogram(hc_col), ColSideColors=clab2, labCol = FALSE,
          #Rowv=NULL, Colv = NULL, ColSideColors = clab3, labCol = FALSE,
          key=TRUE, keysize=1.5, labRow=all_mrs, cexRow=0.35, col=colorspace::diverge_hsv(100),RowSideColors=rlab,
          ColSideColorsSize=3, RowSideColorsSize=1, KeyValueName="Activity Value")
legend("topright",legend=c(unique(cancer_info_vector),"","",unique(icr_type_info_vector),"","",unique(icr_pheno_info_vector)),
       fill=c(unique(Cancers),"black","black",unique(ICR_EorD),"black","black",unique(ICR_Clusters)), border=FALSE, bty="n", y.intersp = 0.7, cex=0.8)
dev.off()



