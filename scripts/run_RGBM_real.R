library(RGBM)
library(doMC)
library(foreach)
library(igraph)
registerDoMC(30)

#setwd('/export/cse/rmall/Network_Analysis/PanCancer_Immunophenotype/Old_Files/ICR_All_Info/scripts/')
setwd('.')

source('get_functions.R')

#Load mechanistic network
load('../Data/Others/me_net_full.Rdata')

#Do the GRN construction and saving
#========================================================================================
filename <- "MESO"
outputpath <- paste0("../Results/",filename)
sample_name <- paste0(filename,"_Full_");
out <- loading_data(filename,M)

V_final <- RGBM(log2(out[[1]]+1), out[[2]], out[[3]], out[[4]], out[[5]], 
                M=5000, nu = 0.001, s_f = 0.3, lf = 1, no_iterations = 3, 
                mink = 0, experimentid = 1, outputpath = outputpath, sample_type = sample_name, real = 1);
write.table(V_final,file=paste("../Results/",filename,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv",sep=""),row.names = T, col.names = T)

# We compress and delete the .csv
system(paste("gzip ../Results/",filename,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv",sep=""))
system(paste("rm ../Results/",filename,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv",sep=""))

