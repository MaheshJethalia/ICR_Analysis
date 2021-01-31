library(corto)
library(doMC)
library(foreach)
library(igraph)
registerDoMC(48)

#setwd('/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Old_Files/ICR_All_Info/scripts/')
setwd('.')

source('get_functions.R')
source('gene-reverse-network.R')

#Load mechanistic network
load('../Data/Others/me_net_full.Rdata')

#Do the GRN construction and saving
#========================================================================================
filename <- "UVM"
outputpath <- paste0("../Results/ARACNE/",filename)
sample_name <- paste0(filename,"_Full_");
out <- loading_data(filename,M)
D <- t(log2(out[[1]]+1))
tfs <- rownames(M)

#Get network using ARACNE-AP algorithm with number of bootstraps to 1000 and perform DPI too
write.table(V_final,file=paste0(outputpath,"/Adjacency_Matrix/Final_Adjacency_Matrix_v2.csv"),row.names = T, col.names = T)

#We compress it and delete the csv
system(paste("gzip ",outputpath,"/Adjacency_Matrix/Final_Adjacency_Matrix_v2.csv",sep=""))
system(paste("rm ",outputpath,"/Adjacency_Matrix/Final_Adjacency_Matrix_v2.csv",sep=""))

