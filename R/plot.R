setwd("/Users/hzr/GitHub/PhyloAcc/R")
source("drawAlign_function.R")
treeData <- prepare_data(tree_path = "../Data/neut_ver3_final.named.mod", species_name = "../Data/species_names.txt", common_name = "../Data/birdname2.txt")

##### 1. read in 
score <- read.table("V2/Combined_elem_lik2.txt", header=T)
score0 <- read.table("Combined_elem_lik_061.txt", header=T)
# order score by loglik_ratio
score <- score[order(-score$log_ratio),]
score$BF2 <- -score$loglik_all + score$loglik_RES
score0 <- score0[order(-score0$log_ratio),]
score0$BF2 <- -score0$loglik_all + score0$loglik_RES  # k=174591,73367
plotZPost(Z, treeData, target_species=NULL, tit=NULL, offset=3)

