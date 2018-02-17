## simulate sequence alignment data under different shift patterns ##
library(MASS)
library(ggplot2)
library(ape)
library(ROCR)

#### if you would like to generated alignments in parallel ####
library(doParallel)
cl<-makeCluster(8)
registerDoParallel(cl)

#### input parameters ####
n=200  # length of subsequence
NE = 1000 # number of elements to simulate
nprior_a = 15  # parameters of Gamma priors of substitution rates
nprior_b = 0.1;
cprior_a = 5
cprior_b = 0.04;
prefix = "../Simulation_ratites/simu_500_200_diffr_2-"  # prefix for the output files



## Read in Q, pi and tree (from *.mod file) ##
alphaB <- c("A","C","G","T")
Q =  matrix(c(-0.891314, 0.108827, 0.630590, 0.151897, 0.150099, -1.149905, 0.130070,0.869737,0.869737,0.130070,-1.149905,
              0.150099, 0.151897, 0.630590, 0.108827, -0.891314),nrow=4, byrow=T)
for(i in 1:4) Q[i,i] <- -sum(Q[i,-i])
eigenvec <- eigen(Q)$vectors
eigenval <-  eigen(Q)$values
eigen_inv <- solve(eigenvec)

pi <- Null(Q);
pi <- pi/sum(pi)

tree <- read.tree("Data/ratite/neut_ver3_final.named.mod") # phylogenetic tree under neutral model
species <- tree$tip.label

nodes_names <- c(tree$tip.label, tree$node.label)

S = length(tree$tip.label)
N = tree$Nnode*2+1
tree$node.label <- 1:(S-1) + S
tree$tip.label <- 1:S


## get rate matrix and adjList ##
branch_mat <- matrix(0, N, N)  # mutation rate from neutral model
adjList <-  rep(list(NULL),N) # Each node, its children in the list 
father <- rep(0,N)

for(i in 1:dim(tree$edge)[1])
{
  branch_mat[tree$edge[i,1], tree$edge[i,2]] <- tree$edge.length[i]
  adjList[[tree$edge[i,1]]] <- c(adjList[[tree$edge[i,1]]], tree$edge[i,2])
  father[tree$edge[i,2]] = tree$edge[i,1]
}

getsubTree <- function(Start = 44){
  
  visited <- rep(0,N) # visited
  S <- Start;  # stack
  
  while(length(S)>0)
  {
    u = S[length(S)];
    S <- S[-length(S)];  # pop
    if(!visited[u])
    {
      visited[u] = 1;
      for(w in adjList[[u]])
      {
        if(visited[w]) next;
        S <- c(S,w);
      }
    }
  }
  return(which(visited==1))
}

#### generate data by DFS ####
## get next base ##
getBase <- function(father,rate, branch=branch_len)
{
  
  P =  eigen_inv * exp(eigenval*rate*branch)
  P <- eigenvec%*%P # rowSums are 1
  P <- P[father,]
  return(sample(1:4,1,prob = P))
}
getData <- function(ind, n, rate, pi, Start = 44){
  if(!is.matrix(rate))
  {
    rate <- matrix(rate,N,N)
  }
  
  rate = rate1[ind]*(1-rate) + rate2[ind]*rate  
  
  # get Data #
  Data <- matrix(0,n,N) # #base * #leaves
  for(i in 1:n)
  {
    tree <- rep(1,N)
    tree[Start] <- sample(1:4,1,prob=pi)  # sample root from pi
    
    visited <- rep(0,N) # visited
    #father <- rep(0,N)
    
    S <- Start;  # stack
    
    while(length(S)>0)
    {
      u = S[length(S)];
      S <- S[-length(S)];  # pop
      #print(u)
      if(!visited[u])
      {
        visited[u] = 1;
        if(u!=Start) tree[u] <- getBase(tree[father[u]],rate[father[u],u], branch_mat[father[u],u])
        for(w in adjList[[u]])
        {
          if(visited[w]) next;
          #father[w] <- u;
          S <- c(S,w);
        }
      }
    }
    Data[i,] <- tree
  }
  for(i in 1:4)
  {
     Data[Data==i] <- alphaB[i]
  }
  return(Data)
}




## different acceleration scenerios for ratites simulations: (2.1) apt branches; (2.2) cas,dro branches; (2.3) rhe branches;  (2.4) strcam; (2.5) all exclude strcam and anoDid; (2.6) all include strcam and anoDid; (2.7) all include 78 subtree; (2.8) randomly selected 5 species to accelerate
nodes_pos <- list(getsubTree(74),getsubTree(76) ,getsubTree(77) ,36,getsubTree(72), c(35, 36,getsubTree(72))) 
nodes_neg <- list(getsubTree(70),sample(c(1:23,31:34), 5))
nodes_all <- c(nodes_pos, nodes_neg)



# generate rate1 and rate2 (i.e. accelerated and conserved rates)
rate1 <- rgamma(NE, nprior_a, scale = nprior_b)
rate2 <- rgamma(NE, cprior_a, scale = cprior_b)
for(i in 1:NE)
{
  if(rate1[i] < rate2[i]) rate1[i] = rate2[i] + 0.5
}

## output bed file (coordinate of each element and rates)##
bedfile <- matrix(0,NE,7)
bedfile[,1] <- 0:(NE-1) # elementID
bedfile[,2] <- seq(0,(NE-1)*n, by = n)
bedfile[,3] <- seq(n,NE*n, by = n)
bedfile[,4] <- 1:NE # elementID
bedfile[,5] = 1; bedfile[,6] = rate2; bedfile[,7] = rate1
options(scipen=10)
write.table(bedfile, file = paste(prefix,1,".bed",sep=""), sep="\t", row.names = F,col.names = F,quote = F)

## output fasta file
for(k in 0:length(nodes_all)) # k=0 generate conserved
{
 
  rate_true <- matrix(1 ,N,N)
  if(k >=1)
  {
    nodes = nodes_all[[k]]
    for(i in nodes )
    {
      j=father[i]
      rate_true[i,j] = rate_true[j,i] = 0 # accelerated branches
    }
  }
  
  registerDoParallel(cl)
  Data <- foreach(i=1:NE, .combine='rbind',.packages="foreach") %dopar% { 
   getData(i, n,rate_true,pi)
  }
  
  ## convert Data to sequence
  colnames(Data) <- nodes_names
  write.dna(t(Data[,1:43]),file = paste(prefix,k,".fasta",sep=""),format = "fasta",nbcol = -1,colsep="")
}





