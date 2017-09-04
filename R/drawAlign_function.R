######## draw acceleration tree and sequence alignment ########
library(seqinr)
library(ggplot2)
library(reshape2)
library(ape)


#### load phylogenetic tree and common name of species ####
## input: tree_path is file of phylogenetic tree; 
#        species_name is a file output by PhyloAcc containing the species name for each column in *postZ* files; 
#        common_name is optional file with three columns: abbreviation of species name appeared in data files and output files; full species names; and species comman name shown on the plot 
## output: a list as input to plotZPost and plotAlign to generate plots
prepare_data <- function(tree_path = "Data/neut_ver3_final.named.mod", species_name = "Data/species_name.txt", common_name = "Data/birdname2.txt")
{
 
    tree <- read.tree(tree_path) # 43 species
    species <- tree$tip.label

    ## reorder species in the output file to be aligned with tree
    species_readin <- read.table(species_name,stringsAsFactors = F)[,1]

    node_idx <- c()
    for(node_name in tree$tip.label) # the last column edge from root ; the first 43 columns are species order the same in ape
    {
      node_idx <-c(node_idx, which(species_readin == node_name))
    }

    for(node_name in tree$node.label) # the last column edge from root ; the first 43 columns are species order the same in ape
    {
      node_idx <-c(node_idx, which(species_readin == node_name))
    }

    ## species common names
    if(!is.null(common_name))
    {
      comnam <- read.table(common_name,sep="\t", row.names = 1, stringsAsFactors = F)
      return(list("tree" = tree, "node_idx" = node_idx, "tip" = comnam[tree$tip.label,2]))
    }
    
    return(list("tree" = tree, "node_idx" = node_idx, "tip" = tree$tip.label))

}

#### Plot the acceleration pattern (posterior of substitution rate on each branch) of one element ####
## input: Z: one row from *postZ* file, posterior of Z for each branch; 0: missing(only for outgroup), 1: neutral, 2: conserved, 3: accelerated
#         tit: title of the figure (BF scores); 
#         treeData: output from  prepare_data
#         target_species: group of species shown in different color (e.g. phenotypically convergent species)
#         offset=2, posterior of Z start from the third column
## output: plot of one element

plotZPost <- function(Z, treeData, target_species=NULL, tit=NULL, offset=3)
{

  species <- treeData$tree$tip.label
  tip.color = rep(1, length(species))
  names(tip.color) <- species
  if(!is.null(target_species)) tip.color[target_species] <- 4
  
  ratio1 = Z[2] # c_rate
  ratio2 = Z[1] # n_rate
  #print(ratio1);print(ratio2);
 

  # color the tree by posterior mean of z
  E = length(Z)
  posterior_z <- rbind( Z[seq(offset + 2,E,by = 4)],Z[seq(offset +3,E,by = 4)], 
                        Z[seq(offset +4,E,by = 4)]) # get posterior of z (only 1???2???3)
  
  # get number of losses by posterior of Z
  #loss = sum(posterior_z[3, 1:length(species)]) - sum(posterior_z[3, -1:-length(species)])
 
  posterior_z <- posterior_z[,treeData$node_idx]
  posterior_z <- posterior_z[,treeData$tree$edge[,2]] # reorder z by edges
  
  mean_z <- colSums(posterior_z * 1:3);  
  
  rbPal <- colorRampPalette(c('gray','purple','springgreen3','firebrick1'))
  mytree= treeData$tree
  mytree$edge.length <-  mytree$edge.length*(posterior_z[1,] + posterior_z[2,] * ratio1 + posterior_z[3,]* ratio2)
  mytree$tip.label <- treeData$tip
  edge_col <- rbPal(91)[round(mean_z*30)+1] 
 
  # color grey is missing
  missing = Z[seq(offset + 1,E,by = 4)]
  tip.color[missing >0] = "azure4"
  plot(mytree,edge.color = edge_col,tip.color = tip.color,edge.width  =3,label.offset = 0.003,no.margin =F, cex=1.1, bg=NA)
  mtext(substitute(paste(t, r[1], "=", r1, ", ", r[2], "=", r2),
    list(t=tit, r1 =round(ratio1,2), r2 = round(ratio2,2))), side = 1, cex=2, line = 0.5)
  
}



#### Plot the sequence alignment for one element ####
## input: kth element to plot
#         align, a matrix of sequence alignment, rows are species, columns are loci
#         bed, a data frame of coordianate of each element
#         legend, “top”(default), “left”,“right”, “bottom”, "none"
## output: heatmap of sequence alignment showing consensus nucleotides, substitutions, indels and N
plotAlign <- function(k, align, bed, treeData, target_species =NULL, legend="top")
{
  
  species <- treeData$tree$tip.label
  cols_sp <- rep(1,length(species))
  names(cols_sp) <- species
  if(!is.null(target_species)) cols_sp[target_species] <- 4
  cols =c("azure2","#df8640","#FFFFFF", "#3794bf") #gainsboro
  
  element1 <- align[,(bed[k+1,2]+1): bed[k+1,3]]
  # get the consensus base pair
  ele_cons <- apply(element1, 2, function(x) { y <- xtabs(~x); 
      cb = names(y)[order(y,decreasing = T)]
      cb <- cb[which(cb%in%c('a','c','g','t'))[1]]
      z <-  rep("subs",length(x)) 
      z[x==cb] <- "cons";
      z[x=='-'] <- "indel";
      z[x=='n'] <- "N"; #| x=='*'
      return(z);
  })
  
  rownames(ele_cons) <- treeData$tip
  dat_m <- melt(ele_cons) 
  p <- ggplot(dat_m, aes(as.factor(Var2), Var1)) + geom_tile(aes(fill = factor(value,levels=c("N","consensus","indel","substitution")))) + 
    scale_fill_manual(values = cols,name="", drop=F) +ylab(NULL) + xlab(paste(ncol(element1), "bp"))+ 
    theme(panel.background = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_text(size=13,colour = cols_sp), 
          axis.ticks.x = element_blank(), axis.text.x =element_blank(),
          axis.title.x = element_text(size = rel(2)),
          axis.title.y = element_text(size = rel(2)),
          legend.text=element_text(size=12),legend.position=legend,
          plot.background = element_rect(fill = "transparent"))
  
  plot(p)
}

