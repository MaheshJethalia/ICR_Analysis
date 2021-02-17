####################################################################
###
### This Script clusters the RNASeq date from TCGA that have been
### normalized and filtered for specific tissue type
### and 1 sample per patient. This script includes a 
### log transformation of the data. Additionally, optimal kalinsky
### is calculated.
### 
### Input data:
### RNAseq data- 
### "../Data/,"Cancer,"/",Cancer,"_gene_RNASeq_normalized_TP_filtered.Rdata"

### Output data are saved as Rdata file:
### "../Data/",Cancer,"/",Cancer,"_ICR_cluster_assignment_k2-6.Rdata"
### which includes: (1) table_cluster_assignment and (2) optimal.calinsky.
#####################################################################

rm(list=ls())

setwd(".")                                                                    # Setwd to location were output files have to be saved.

source("ipak.function.R")
source("get_functions.R")

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper","ConsensusClusterPlus")
required.bioconductor.packages = c("clue")
ipak(required.packages)
ibiopak(required.bioconductor.packages)


# Set Parameters
Cancer = "BLCA"               
# Load ICR genes data
load("../Data/Others/ICR_genes.RData") 

#Load RNA-Seq data
out <- load(paste0("../Data/",Cancer,"/",Cancer,"_gene_RNAseq_normalized_TP_filtered.Rdata"))
D <- t(filtered.norm.RNAseqData)
ICR_subset_RNAseq = D[,colnames(D) %in% ICR_genes]
ICR_subset_RNAseq_log2 = log(ICR_subset_RNAseq +1, 2) 


#Perform the consensus clustering
ddist = dist(ICR_subset_RNAseq_log2)

ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                               maxK = 6,                                                                              # set K
                                               pItem = 0.8,
                                               reps=5000,                                                                             # set repeats
                                               title=paste0(Cancer,".ICR.reps5000"),                    # Output filename
                                               clusterAlg = "hc",                                                                     # clustering Algorithm : Hierarchiocal clustering
                                               innerLinkage = "ward.D2",                                                              # for color coding the clusters use tmyPal = ...
                                               finalLinkage = "complete",
                                               plot = 'pdf',                                                                          # write resut to pdf (Alt:png)
                                               writeTable = TRUE,
                                               verbose = TRUE)
outputfiles = list.files(paste0(Cancer,".ICR.reps5000"), full.names = TRUE)
class_files = outputfiles[grep("consensusClass", outputfiles)]

N.files = length(class_files)
table_cluster_assignment = data.frame(ICRscore = rowMeans(ICR_subset_RNAseq_log2))

for (j in 1:N.files){
  file = paste0("./", class_files[j])
  consensus_class = read.csv(file = file,header=FALSE)
  group = paste0("Group_k",j+1)
  colnames(consensus_class) = c("PatientID", group)
  rownames(consensus_class) = consensus_class$PatientID
  consensus_class$PatientID = NULL
  table_cluster_assignment[,group] = consensus_class[,group][match(rownames(table_cluster_assignment), rownames(consensus_class))]
  
  transl_table_ICR_cluster = aggregate(ICRscore~get(group),data = table_cluster_assignment, FUN=mean)
  colnames(transl_table_ICR_cluster) = c(group,"mean_ICRscore")
  transl_table_ICR_cluster = cbind(transl_table_ICR_cluster[order(transl_table_ICR_cluster$mean_ICRscore),],ICR_name=paste0("ICR",c(1:(j+1))))
  
  ICR_cluster = paste0("ICR_cluster_k",j+1)
  table_cluster_assignment[,ICR_cluster] = transl_table_ICR_cluster$ICR_name[match(table_cluster_assignment[,group],
                                                                                   transl_table_ICR_cluster[,group])]


  
  
  #calinsky
  sHc <- hclust(ddist, method = "ward.D2")
  aCalinsky <- calinsky(sHc,gMax=10)
  # pdf(file = paste0("Data/", Cancer, "/", Cancer, "_ICR_cluster_assignment_k2-6.Calinsky.pdf"), width = 16, height = 6)
  # plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
  # text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
  # dev.off()
  optimal.calinsky = which(aCalinsky == max(aCalinsky[3:5]))
  
  #save data
  outputname = paste0("../Data/", Cancer, "/", Cancer, "_ICR_cluster_assignment_k2-6_v2.Rdata")
  save(table_cluster_assignment,optimal.calinsky, file = outputname)
}
system(paste0("rm -rf ",Cancer,"*"))
