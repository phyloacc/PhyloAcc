# PhyloAcc
This is a software to detect the shift pattern of DNA substitution rate of a genomic region, and identify genomic elements accelerated in some specific species from a set of conserved elements. The underlying model assumes a latent discrete state (Z) of relative substitution rate of each branch on the phylogeny which can be neutral, conserved and accelerated. For each genomic element, it will start from neutral rate at the common ancestor of the phylogeny, transit to conserved state and then reach accelerated state in some lineages. Our method utilizes adaptive collapsed Gibbs sampling to obtain the pattern of substitution rate shifts (posterior distribution of Z) as well as relative substitution rates of conserved and accelerated state. To identify DNA elements with accelerating on specific branches, we compared marginal likelihood under three models: null model (M0) where all branches are neutral or conserved; accelerated model (M1) in which branches leading to target species are accelerated; and full model (M2), with no constraints on latent states. Then we used two Bayes factors: between M1 and M0 (BF1) and between M1 and M2 (BF2) as criteria to identify DNA elements accelerated exclusively in target lineages.

We need some prelimenary inputs which might be generated from other softwares: (1) a Phylogeny in mod format. The file can be output from phyloFit in PHAST package, with the transition rate matrix of bases and branch lengths. Our model assumes that those are the substitution rates under neutral evoluation and will estimate conserved and accelerated rate relatived to the neutral rate. In our study, we used genome-wide four-fold degenerate sites to estimtate the rate matrix and branch lengths. (2) a multiple alignment file concatenating sequences of all input conserved elements in FASTA format, and (3) a bed file with the position of each individual element in the coordinate of concatenated alignment file (0-based).

We also need a parameter file, which contains the pathes for input files and output directory, information of species and parameters for MCMC. In the parameter file, each parameter is specified in a line with the parameter name at the beginning followed by paremeter value. The parameters are: 
PHYTREE_FILE: the path of Phylogeny (.mod)  
SEG_FILE: the path of bed file for genomic regions
ALIGN_FILE: the path of multiple alignment file (.fasta)
RESULT_PREFIX: the prefix for output files
BURNIN: number of initial iterations to discard before equilibrium of the chain (default: 500)
MCMC: number of MCMC iterations (default: )
CHAIN: 1
TARGETSPECIES: species of interest. E.g. species potentially lost conservation or with convergent phenotype changes.
OUTGROUP: outgroup species of the phylogeny. These species are not considered to be accelerated in our model. 
CONSERVE: species assumed to be mostly conserved. The algorithm will filter out elements with alignment gaps in more than 50% of the conserved species. 
REF: galGal
INDEL 0
GAPCHAR -
