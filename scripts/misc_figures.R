library(data.table)
library(ggplot2)
library(gplots)
library(Matrix)
library(erer)
library(doMC)
library(heatmap3)
library(mgcv)
library(devtools)
library(gprofiler2)


#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

setwd('/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Old_Files/ICR_All_Info/scripts/')
setwd('.')

source('gene-reverse-network.R')
source('get_functions.R')

convert_character_vector <- function(df)
{
  for (i in 1:ncol(df))
  {
    df[,i] <- as.character(as.vector(df[,i]))
  }
  return(df)
}

#For survival analysis
###################################################################################################################
all_cold_mrs_df <- read.table("../Results/Revised_Text_Results/All_MRS_ICR_Low_Supplementary_Table_S7d.csv",sep="&",header=TRUE)
all_cold_mrs_df <- all_cold_mrs_df[all_cold_mrs_df$Median_ICR_Enabled>=0 & all_cold_mrs_df$Median_ICR_Disabled<0,]
all_cold_mrs <- as.character(as.vector(all_cold_mrs_df[,1]))

all_hot_mrs_df <- read.table("../Results/Revised_Text_Results/All_MRS_ICR_High_Supplementary_Table_S7c.csv",sep="&",header=TRUE)
all_hot_mrs_df <- all_hot_mrs_df[all_hot_mrs_df$Median_ICR_Enabled<0 & all_hot_mrs_df$Median_ICR_Disabled>=0,]
all_hot_mrs <- as.character(as.vector(all_hot_mrs_df[,1]))
all_mrs <- union(all_hot_mrs,all_cold_mrs)

RGBM_path <- c("../Data/");

icr_disabled_cancers <- c("LGG","PAAD","UVM","KIRC")
icr_enabled_cancers <- c("BLCA","BRCA","HNSC","LIHC","SARC","STAD","SKCM","UCEC")

cancers <- c(icr_disabled_cancers,icr_enabled_cancers)
for (i in 1:length(cancers))
{
  cancer_type <- cancers[i]
  load(paste0("../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
  
  high_indices_table <- read.table(paste0("../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_high_indices.csv"),header=TRUE)
  low_indices_table <- read.table(paste0("../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_low_indices.csv"),header=TRUE)
  amat_to_use <- amat[,c(high_indices_table$x,low_indices_table$x)]
  amat_to_use[amat_to_use>0] <- amat_to_use[amat_to_use>0]/max(amat_to_use)
  amat_to_use[amat_to_use<0] <- amat_to_use[amat_to_use<0]/abs(min(amat_to_use))
  samples <- dim(amat_to_use)[2]
  high_index_last <- nrow(high_indices_table)
  low_index_last <- ncol(amat_to_use)
  valid_mrs <- rownames(amat_to_use)[rownames(amat_to_use) %in% all_mrs]
  icr_pheno_info_vector <- c(rep("ICR High",nrow(high_indices_table)),
                             rep("ICR Low",nrow(low_indices_table)))
  sample_ids <- colnames(amat_to_use)
  valid_amat_to_use <- amat_to_use[valid_mrs,]
  save(list = c('valid_amat_to_use','sample_ids','icr_pheno_info_vector'), file=paste0("../Results/PanCancer/",cancer_type,"_activity_info_for_survival.Rdata"))
}

#Analyze the PRECOG dataset
##########################################################################################################
#predictions_df <- readRDS("/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Data/PRECOG/predictors_all_merged3.rds")
#expr_list1 <- readRDS("/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Data/PRECOG/es_list1.rds")
#expr_list2 <- readRDS("/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Data/PRECOG/es_list2.rds")
#mapping_df <- read.table("../Data/PRECOG/Mapping_GEO.csv",header=TRUE,sep=",")
#mapping_df <- convert_character_vector(mapping_df)

# #Get only those samples for which cancer type is available
# cancer_list <- NULL
# rev_predictions_df <- predictions_df[predictions_df$Study %in% mapping_df$Study,]
# for (i in 1:nrow(rev_predictions_df))
# {
#   study_type <- rev_predictions_df$Study[i]
#   cancer_type <- mapping_df[mapping_df$Study==study_type,]$Cancer  
#   cancer_list <- c(cancer_list,cancer_type)
# }
# rev_predictions_df$Cancer <- cancer_list
# 
# #Get largest cohort for each cancer type
# geo_cancer_mat <- as.matrix(table(rev_predictions_df$Study,rev_predictions_df$Cancer))
# unique_cancers <- colnames(geo_cancer_mat)
# geo_cancer_edgelist <- NULL
# for (i in 1:length(unique_cancers))
# {
#   id <- which.max(geo_cancer_mat[,colnames(geo_cancer_mat)==unique_cancers[i]])
#   geo_id <- rownames(geo_cancer_mat)[id]
#   sample_size <- geo_cancer_mat[geo_id,unique_cancers[i]]
#   temp <- cbind(geo_id,unique_cancers[i],sample_size)
#   geo_cancer_edgelist <- rbind(geo_cancer_edgelist,temp)
# }
# geo_cancer_edgelist <- as.data.frame(geo_cancer_edgelist)
# colnames(geo_cancer_edgelist) <- c("GEO_Accession_Id","Cancer","Sample_Size")
# geo_cancer_edgelist$GEO_Accession_Id <- as.character(as.vector(geo_cancer_edgelist$GEO_Accession_Id))
# geo_cancer_edgelist$Cancer <- as.character(as.vector(geo_cancer_edgelist$Cancer))
# geo_cancer_edgelist$Sample_Size <- as.numeric(as.vector(geo_cancer_edgelist$Sample_Size))
# 
# geo_cancer_edgelist <- geo_cancer_edgelist[geo_cancer_edgelist$Sample_Size>100,]
# 
# final_expr_list <- list()
# k <- 1
# for (i in 1:nrow(geo_cancer_edgelist))
# {
#   geo_id <- geo_cancer_edgelist$GEO_Accession_Id[i]
#   if (sum(names(expr_list1) %in% geo_id)>0)
#   {
#     id <- which(names(expr_list1)==geo_id)
#     final_expr_list[[k]] <- expr_list1[[id]]
#   }
#   else if (sum(names(expr_list2) %in% geo_id)>0)
#   {
#     id <- which(names(expr_list2)==geo_id)
#     final_expr_list[[k]] <- expr_list2[[id]]
#   }
#   k <- k + 1
# }
# names(final_expr_list) <- geo_cancer_edgelist$GEO_Accession_Id
# 
# #Create the expression dataset and ICR High and Low indices per cancer type
# for (i in 1:length(final_expr_list))
# {
#   #Get GEO Id and Cancer Id
#   geo_id <- names(final_expr_list[i])
#   cancer_type <- geo_cancer_edgelist[geo_cancer_edgelist$GEO_Accession_Id==geo_id,]$Cancer
#   
#   #Get ICR Scores and sample ids of samples in dataset
#   icr_scores <- rev_predictions_df[rev_predictions_df$Study==geo_id & rev_predictions_df$Cancer==cancer_type,]$ICR
#   sample_ids <- as.character(rev_predictions_df[rev_predictions_df$Study==geo_id & rev_predictions_df$Cancer==cancer_type,]$ID)
#   
#   #Get ICR High and ICR Low cutoffs
#   quantiles_icr_scores <- quantile(icr_scores)
#   icr_low_cutoff <- quantiles_icr_scores[[2]]
#   icr_high_cutoff <- quantiles_icr_scores[[4]]
#   
#   #Get ICR High and ICR Low scores
#   icr_low_ids <- which(icr_scores<icr_low_cutoff)
#   icr_high_ids <- which(icr_scores>icr_high_cutoff)
#   
#   #Get the expression matrix from GEO accession
#   if (i==3)
#   {
#     sample_names <- strsplit(colnames(final_expr_list[[i]]),"_")
#     rev_sample_names <- NULL
#     for (l in 1:length(sample_names))
#     {
#       rev_sample_names <- c(rev_sample_names,sample_names[[l]][1])
#     }
#     colnames(final_expr_list[[i]]) <- rev_sample_names
#   }
#   expr_matrix <- exprs(final_expr_list[[i]])[,sample_ids]
#   
#   #Map gene ids to gene names (HGNC symbols)
#   gene_ids <- rownames(expr_matrix)
#   
#   if (i==6)
#   {
#     rev_gene_ids <- rep(0,length(gene_ids))
#     for (l in 1:length(gene_ids))
#     {
#       gene_id <- unlist(strsplit(gene_ids[l]," "))[1]
#       rev_gene_ids[l] <- as.numeric(gene_id)
#     }
#     query_out <- gconvert(query= rev_gene_ids, organism = "hsapiens", numeric_ns="ENTREZGENE_ACC",
#                           target = "HGNC", mthreshold = Inf, filter_na = TRUE)
#     gene_ids <- rev_gene_ids
#   }else
#   {
#     query_out <- gconvert(query = gene_ids, organism = "hsapiens",
#                           target="HGNC", mthreshold = Inf, filter_na = TRUE)
#   }
#   
#   gene_names <- NULL
#   for (j in 1:length(gene_ids))
#   {
#     if (gene_ids[j] %in% query_out$input)
#     {
#       gene_name <- query_out[query_out$input==gene_ids[j],]$name[1]
#     }
#     else{
#       gene_name <- "NA"
#     }
#     gene_names <- c(gene_names,gene_name)
#   }
#   gene_indices <- which(gene_names!="NA")
#   rev_expr_matrix <- expr_matrix[gene_indices,]
#   rownames(rev_expr_matrix) <- gene_names[gene_indices]
#   
#   #Get revised expression matrix with ordered ICR High and ICR Low
#   D <- rev_expr_matrix[,c(icr_high_ids,icr_low_ids)]
#   load(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_Correlation_matrix.Rdata"))
#   load(paste0("Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_regulon_FGSEA.Rdata"))
#   
#   #Revise the regulons with TFs present in expression matrix and targets present in expression matrix
#   all_tfs <- names(regulon_sets)
#   all_targets <- rownames(D)
#   rev_tfs <- all_tfs[all_tfs %in% all_targets]
#   rev_regulon_sets <- list()
#   k <- 1
#   for (tf in rev_tfs)
#   {
#     id <- which(all_tfs==tf)
#     rev_targets <- regulon_sets[[id]][regulon_sets[[id]] %in% all_targets]
#     rev_regulon_sets[[k]] <- rev_targets
#     k <- k+1
#   }
#   names(rev_regulon_sets) <- rev_tfs
#   rev_corr_matrix <- corr_matrix[,which(colnames(corr_matrix) %in% all_targets)]
#   
#   #Generate the activity matrix
#   amat <- activity_mc(mexp = D, cormat = rev_corr_matrix, tflist = names(rev_regulon_sets), tau = 0.0)
#   save(list=c('amat'),file=paste0("Data/PRECOG/",cancer_type,"/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
#   icr_info <- c(rep("ICR High",length(icr_high_ids)),rep("ICR Low",length(icr_low_ids)))
#   write.table(icr_info,paste0("Data/PRECOG/",cancer_type,"/",cancer_type,"_ICR_Info.csv"),row.names=T,col.names=F,quote=F,sep=",")
# }

#Get the activity from multiple PRECOG datasets
########################################################################################################################
precog_cancers <- c("BLCA","BRCA","COAD","GBM","HNSC","LUAD","OV","SKCM")
icr_type <- c("ICR Enabled","ICR Enabled","ICR Neutral","ICR Neutral","ICR Enabled","ICR Neutral","ICR Neutral","ICR Enabled")
inputpath <- "../Data/PRECOG/"

#Get all MRs of interest
all_cold_mrs_df <- read.table("../Results/Revised_Text_Results/All_MRS_ICR_Low_Supplementary_Table_S7d.csv",sep="&",header=TRUE)
all_cold_mrs_df <- all_cold_mrs_df[all_cold_mrs_df$Median_ICR_Enabled>0 & all_cold_mrs_df$Median_ICR_Disabled>0,]
all_cold_mrs <- as.character(as.vector(all_cold_mrs_df[,1]))

all_hot_mrs_df <- read.table("../Results/Revised_Text_Results/All_MRS_ICR_High_Supplementary_Table_S7c.csv",sep="&",header=TRUE)
all_hot_mrs_df <- all_hot_mrs_df[all_hot_mrs_df$Median_ICR_Enabled>=0 & all_hot_mrs_df$Median_ICR_Disabled>=0,]
all_hot_mrs <- as.character(as.vector(all_hot_mrs_df[,1]))
all_mrs <- c(all_hot_mrs,all_cold_mrs)

cancer_info_vector <- NULL
icr_pheno_info_vector <- NULL
icr_type_info_vector <- NULL
sample_ids <- NULL

#Get total no of cancer samples
cancer_samples <- 0
for (i in 1:length(precog_cancers))
{
  cancer_type <- precog_cancers[i]
  table_indices <- read.table(paste0(inputpath,cancer_type,"/",cancer_type,"_ICR_Info.csv"),header=FALSE,sep=",")
  cancer_samples <- cancer_samples+nrow(table_indices)
}

activity_matrix <- Matrix(0,nrow=length(all_mrs),ncol=cancer_samples)
rownames(activity_matrix) <- all_mrs
counter=0;

#Put values in the activity matrix
for (i in 1:length(precog_cancers))
{
  cancer_type <- precog_cancers[i]
  load(paste0(inputpath,cancer_type,"/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
  table_indices <- read.table(paste0(inputpath,cancer_type,"/",cancer_type,"_ICR_Info.csv"),header=FALSE,sep=",")
  table_indices$V2 <- as.character(as.vector(table_indices$V2))
  high_index_last <- sum(table_indices$V2=="ICR High")
  low_index_last <- nrow(table_indices)
    
  amat_to_use <- amat  
  amat_to_use[amat_to_use>0] <- amat_to_use[amat_to_use>0]/max(amat_to_use)
  amat_to_use[amat_to_use<0] <- amat_to_use[amat_to_use<0]/abs(min(amat_to_use))
  samples <- dim(amat_to_use)[2]
  
  cancer_info_vector <- c(cancer_info_vector,rep(cancer_type,samples))
  icr_pheno_info_vector <- c(icr_pheno_info_vector,rep("ICR High",sum(table_indices$V2=="ICR High")),
                             rep("ICR Low",sum(table_indices$V2=="ICR Low")))
  icr_type_info_vector <- c(icr_type_info_vector,rep(icr_type[i],samples))
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

colnames(activity_matrix) <- sample_ids

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
activity_matrix <- as.matrix(activity_matrix)
#2nd column is ICR High vs ICR Low
hc.cols1 <- order(clab[,2],decreasing = T)
clab1 <- clab[hc.cols1,]
activity_matrix <- activity_matrix[,hc.cols1]
newhc.cols <- NULL
#3rd column is ICR Enabled, ICR Disabled and ICR Neutral clusters
for (i in 1:length(unique(clab[,2])))
{
  color_id <- unique(clab[,2])[i]
  newhc.cols <- c(newhc.cols,length(newhc.cols)+order(clab1[clab1[,2]==color_id,3],decreasing = T))
}
activity_matrix <- activity_matrix[,newhc.cols]
clab2 <- clab1[newhc.cols,]
lower_val <- quantile(activity_matrix,probs=c(0.05))
upper_val <- quantile(activity_matrix,probs=c(0.95))
activity_matrix[activity_matrix<lower_val] <- lower_val
activity_matrix[activity_matrix>upper_val] <- upper_val

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
hc.cols.icr_high <- hclust(dist(t(activity_matrix[,all_icr_high_indices])),method='ward.D2')
hc.cols.icr_low <- hclust(dist(t(activity_matrix[,all_icr_low_indices])),method='ward.D2')
hc.cols <- c(hc.cols.icr_high$order,(length(all_icr_high_indices)+hc.cols.icr_low$order))
hc_col_list <- list(hc.cols.icr_high,hc.cols.icr_low)
hc_col <- .merge_hclust(hc_col_list)
#clab3 <- clab2[hc.cols,]
#Activity_Matrix <- Activity_Matrix[,hc.cols]

all_hot_mr_indices <- which(rlab[1,]==icr_colors[1])
all_cold_mr_indices <- which(rlab[1,]==icr_colors[2])
hc.rows.icr_high <- hclust(dist(activity_matrix[all_hot_mr_indices,]),method='ward.D2')
hc.rows.icr_low <- hclust(dist(activity_matrix[all_cold_mr_indices,]),method='ward.D2')
hc.rows <- c(hc.rows.icr_high$order,(length(all_hot_mr_indices)+hc.rows.icr_low$order))
hc_row_list <- list(hc.rows.icr_high,hc.rows.icr_low)
hc_row <- .merge_hclust(hc_row_list)

#Main Plotting Function
pdf("../Results/Revised_Figures/Pdfs/Heatmap_All_MRs_PRECOG_Figure_Revised(v2)_5C.pdf",height = 13, width=15, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0)
p1 <- heatmap.3(activity_matrix, scale="none", dendrogram="both", margins=c(4,20),
          Rowv=as.dendrogram(hc_row), Colv=as.dendrogram(hc_col), ColSideColors=clab2, labCol = FALSE,
          #Rowv=NULL, Colv = NULL, ColSideColors = clab3, labCol = FALSE,
          key=TRUE, keysize=1.5, labRow=all_mrs, cexRow=0.35, col=colorspace::diverge_hsv(100),RowSideColors=rlab,
          ColSideColorsSize=3, RowSideColorsSize=1, KeyValueName="Activity Value")
legend("topright",legend=c(unique(cancer_info_vector),"","",unique(icr_type_info_vector),"","",unique(icr_pheno_info_vector)),
       fill=c(unique(Cancers),"white","white",unique(ICR_EorD),"white","white",unique(ICR_Clusters)), border=FALSE, bty="n", y.intersp = 0.7, cex=0.8)
dev.off()

ordered_mrs <- rev(rownames(activity_matrix)[p1$rowInd])
high_ids <- clab2[,2]=="yellow"
low_ids <- clab2[,2]=="green"
activity_df <- perform_wilcox_test(activity_matrix[ordered_mrs,high_ids],activity_matrix[ordered_mrs,low_ids])
activity_df <- activity_df[,c(2:ncol(activity_df))]
activity_df$Mean1 <- round(activity_df$Mean1,3)
activity_df$Mean2 <- round(activity_df$Mean2,3)
activity_df$FC_Mean <- round(activity_df$FC_Mean,3)
write.table(activity_df,"../Results/Revised_Text_Results/Diff_MRS_Activity_PRECOG_Revised(v2)_Table_13.csv",row.names=T,col.names=T,quote=F,sep=",")
