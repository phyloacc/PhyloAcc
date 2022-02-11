## simulate sequence alignment data under different shift patterns and gBGC##
library(MASS)
library(ggplot2)
library(ape)
library(ROCR)

#### if you would like to generated alignments in parallel ####
library(doParallel)
cl<-makeCluster(8)
registerDoParallel(cl)

## input parameter ##
n=1000  # length of sequence
NE = 200 # number of elements

nprior_a = 6 # parameters for gamma prior of selection coefficients and gBGC coefficients
nprior_b = 0.4;
cprior_a = 5
cprior_b = 0.4;


prefix = "../V2_GBGC/Simulation/simu_200_200_diffr_2-" #prefix for the output files

#### Read in Q, pi and tree ####
alphaB <- c("A","C","G","T")
Q =  matrix(c(-0.891314, 0.108827, 0.630590, 0.151897, 0.150099, -1.149905, 0.130070,0.869737,0.869737,0.130070,-1.149905,
              0.150099, 0.151897, 0.630590, 0.108827, -0.891314),nrow=4, byrow=T)
for(i in 1:4) Q[i,i] <- -sum(Q[i,-i])
eigenvec <- eigen(Q)$vectors
eigenval <-  eigen(Q)$values
eigen_inv <- solve(eigenvec)

B = matrix(0, 4, 4)  # biased Gene conversion matrix
B[c(1,1,4,4), c(2,3,2,3)] = 1
B[c(2,3,2,3), c(1,1,4,4)] = -1


pi <- Null(Q);
pi <- pi/sum(pi)

tree <- read.tree("Data/ratite/neut_ver3_final.named.mod") # 43 species
species <- tree$tip.label
nodes_names <- c(tree$tip.label, tree$node.label)

S = length(tree$tip.label)
N = tree$Nnode*2+1
tree$node.label <- 1:(S-1) + S
tree$tip.label <- 1:S

## get rate matrix and adjList##
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
  #father <- rep(0,N)
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
getEigen <- function(s, b)
{
  Y = (s + b*B)/(1 - exp(-(s + b * B)))
  Y[is.na(Y)] = 1
  P = Q * Y
  for(i in 1:4) P[i,i] <- -sum(P[i,-i])
  eigenvec <- eigen(P)$vectors
  eigenval <-  eigen(P)$values
  eigen_inv <- solve(eigenvec)
  return(list(eigenvec,eigenval,eigen_inv ))
}

## get next base ##
getBase <- function(father,decomp, b=0, rate=0, branch)
{
  if(b==0)
  {
    P =  eigen_inv * exp(eigenval*rate *branch)
    P <- eigenvec%*%P # rowSums are 1
  }else{
    
    P =  decomp[[3]] * exp(decomp[[2]]*branch)
    P <- decomp[[1]]%*%P # rowSums are 1
  }
  
  P <- P[father,]
  return(sample(1:4,1,prob = P))
}

#ind: index of simulated elements; n: length of sequence; rateInd: binary matrix for shift pattern; bias_node: vector indicating gBGC effect for each branch
getData <- function(ind, n, rateInd, bias_node, pi, Start = 44){
  s = rate1[ind]*(1-rateInd) + rate2[ind]*rateInd  # get rate from proportion
  rate = s/(1 - exp(-s)) # rate with B=0
  
  # do eigen decomp of rate matrix
  decomp <- list()
  decomp[[1]] = getEigen(rate1[ind], bias[ind])
  decomp[[2]] = getEigen(rate2[ind], bias[ind])
  
  
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
      
      if(!visited[u])
      {
        visited[u] = 1;
        
        if(u!=Start){
          sInd = rateInd[father[u],u] + 1
          tree[u] <- getBase(tree[father[u]],decomp[[sInd]],bias_nodes[u],
                             rate[father[u],u], branch_mat[father[u],u])
        }
        for(w in adjList[[u]])
        {
          if(visited[w]) next;
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


## different acceleration scenerios for ratites simulations:
##(2.1) apt branches; (2.2) cas,dro branches; (2.3) rhe branches;  (2.4) strcam; (2.5) all exclude strcam and anoDidl (2.6) all include strcam and anoDid; (2.7) all include 78 subtree; (2.8) randomly selected 5 species to accelerate
nodes_pos <- list(getsubTree(74),getsubTree(76) ,getsubTree(77) ,36,getsubTree(72), c(35, 36,getsubTree(72)))   
nodes_neg <- list(getsubTree(70),sample(c(1:23,31:34), 5))
nodes_all <- c(nodes_pos, nodes_neg)


# generate rate1 and rate2
rate1 <- rgamma(NE, nprior_a, scale = nprior_b) -2 
rate2 <- -rgamma(NE, cprior_a, scale = cprior_b)

# generate gBGC coefficients Gamma(7, 0.5)
bias <-  rgamma(NE, 7, scale = 0.5)

rate1 <- rep(0, NE) 
rate2 <- rep(-3, NE)
bias <-  runif(NE, 0, 2) #all ratites acceleration, s: -2 vs. 0; 3-0: low gBGC, 3-1:median gBGC; 3-2:high gBGC 

bias_nodes <- rep(1, N)

## different gBGC scenerios:
#nodes_bias <- c(c(27,28),sample(c(1:23,31:34), 3))   # getsubTree(77) and 3 randomly branches
bias_nodes <- rep(0, N)
bias_nodes[nodes_bias] = 1

#rate1[rate1<1] = 1
for(i in 1:NE)
{
  if(rate1[i] < rate2[i]) rate1[i] = rate2[i] + 1
}
bedfile <- matrix(0,NE,8)
bedfile[,1] <- 0:(NE-1) # elementID
bedfile[,2] <- seq(0,(NE-1)*n, by = n)
bedfile[,3] <- seq(n,NE*n, by = n)
bedfile[,4] <- 1:NE # elementID
bedfile[,5] = bias; bedfile[,6] = rate2; bedfile[,7] = rate1
bedfile[,8] = paste(nodes_bias, collapse = ",") #
options(scipen=10)
write.table(bedfile, file = paste(prefix,1,".bed",sep=""), sep="\t", row.names = F,col.names = F,quote = F)


for(k in c(0:length(nodes_pos))) # k=0 generate conserved, all
{
  rate_true <- matrix(1 ,N,N)
  if(k >=1)
  {
    nodes = nodes_pos[[k]] #nodes_bias
    for(i in nodes )
    {
      j=father[i]
      rate_true[i,j] = rate_true[j,i] = 0
    }
  }

  registerDoParallel(cl)
  Data <- foreach(i=1:NE, .combine='rbind',.packages="foreach") %dopar% { 
   getData(i, n,rate_true,bias_nodes, pi)
  }
  
  ## convert Data to sequence
  colnames(Data) <- nodes_names
  write.dna(t(Data[,1:43]),file = paste(prefix,2,".fasta",sep=""),format = "fasta",nbcol = -1,colsep="")
}






