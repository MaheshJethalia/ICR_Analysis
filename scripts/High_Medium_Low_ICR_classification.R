#################################################################
###
### This script re-classifies patients to High-Medium-Low ICR
### clusters as specified in 

### Input- & Outputfiles:
### "../Data/",Cancer,"/",Cancer, "_ICR_cluster_assignment_k2-6.Rdata"
### (Rdata files get updated using this script)
#################################################################

# Setup environment
rm(list=ls())

setwd(".")                                                                    # Setwd to location were output files have to be saved.
source("ipak.function.R")

required.packages = c("plyr")
ipak(required.packages)

# Set Parameters
Cancer = "BLCA"                                                                                                     # Specify the cancertypes

# Load data
Cluster_file = paste0("../Data/",Cancer,"/",Cancer,"_ICR_cluster_assignment_k2-6_v2.Rdata")
load(Cluster_file)

if(optimal.calinsky == 3){
  table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$ICR_cluster_k3)
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR1"] = "ICR Low"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR2"] = "ICR Medium"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR High"
}
if(optimal.calinsky == 4){
  table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$ICR_cluster_k4)
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR1"] = "ICR Low"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR2"] = "ICR Medium"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR Medium"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR4"] = "ICR High"
}
if(optimal.calinsky == 5){
  table_cluster_assignment$HML_cluster = as.character(table_cluster_assignment$ICR_cluster_k5)
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR1"] = "ICR Low"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR2"] = "ICR Medium"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR3"] = "ICR Medium"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR4"] = "ICR Medium"
  table_cluster_assignment$HML_cluster[table_cluster_assignment$HML_cluster == "ICR5"] = "ICR High"
}

save(table_cluster_assignment,optimal.calinsky, file = Cluster_file)
