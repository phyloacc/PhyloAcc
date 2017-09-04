source("drawAlign_function.R")
treeData <- prepare_data(tree_path = "../Data/neut_ver3_final.named.mod", species_name = "../Data/species_name.txt", common_name = "../Data/birdname.txt")
plotZPost <- function(Z, treeData, target_species=NULL, tit=NULL, offset=3)

