# Parameter file
In the parameter file, each parameter is specified in a line with the parameter name at the beginning followed by the parameter value. The parameters are: 
* **Input and output**: 
  * *PHYTREE_FILE*: the path of phylogeny (.mod)  
  * *SEG_FILE*: the path of bed file for genomic regions. At least 3 columns in the bed file. You could have more columns in this file but the program will only read in the first 3 columns. The first column is element name or ID (different from usual bed file), the second and third columns are the starting and ending positions of each element (in the coordinate of the whole multiple alignment). The program assumes that the alignment file concatenates all the elements together and will only use the second and third columns in the bed file. If concatenating multiple chromosomes, the coordinate of elements on the current chromosome should not start from zero but should add to the previous chromosome. The program will internally generated a No. for each element which is the order in the input bed file, and it will use No. in the outputs and plot functions.
  * *ALIGN_FILE*: the path of multiple alignment file (.fasta). The name of the species in the alignment file has to the same as in the phylogenetic tree!
  * *RESULT_FOLDER*: the output folder. The folder should exist.
  * *PREFIX*: the prefix for output files (default: test).
  * *ID_FILE* (optional): only compute elements in this file. (The element is numbered by its order in the input bed file starting from 0). If not specified, the program will compute all elements in the input file.  
  * *VERBOSE*: 0 or 1. If it's 1, the algorithm will output some intermediate results to console and MCMC trace for each element (default: 0). Should set to 0 if computing many elements, otherwise the output file is too large. 
 
* **Specify species on the phylogeny**:
  * *TARGETSPECIES*: species of interest. E.g. species potentially lost conservation or with convergent phenotype changes.
  * *OUTGROUP*: outgroup species of the phylogeny. These species are not considered to be accelerated in our model. 
  * *CONSERVE*: species assumed to be mostly conserved. The algorithm will filter out elements "missing" in more than *CONSERVE_PROP* of the conserved species. Input conserved species should exclude target species.
  * *CONSERVE_PROP*: filter out elements "missing" in more than *CONSERVE_PROP* of the conserved species (default: 0.8).
  * *PRUNE_TREE*: Whether to prune "missing" branches besides outgroup. The "missing" branches inside outgroup are always pruned. (default: false, no prune).
  
* **Alignment Gaps and Filtering**:  
  * *GAPCHAR*: the character for alignment gaps. (default: -). Should be one char.
  * *GAP_PROP*: if the sequence alignment of a species contains gaps for more than *GAP_PROP* of the whole element, then we say that the element is "missing" in that species (default: 0.8).   
  * *TRIM_GAP_PERCENT*: Trim the loci with indels or unknown base pairs in more than *TRIM_GAP_PERCENT* of all species. (default: 1, no trim).
  * *CONSTOMIS*: the probability of "missing" under the conserved state. Should be small (default: 0.01). 
  * *MIN_LEN*: The trimmed element with length less than *MIN_LEN* will be filtered out. (default: 50).

* **(Hyper)Parameters and initial values**:
  * *INIT_GRATE*: the initial transition probability from background to conserved state (default: 0.5).
  * *INIT_LRATE*: the initial transition probability from conserved to accelerated state (default: 0.3).
  * *HYPER_LRATE_A, HYPER_LRATE_B*: the parameters for the beta prior of loss probability (default: 1,9).
  * *HYPER_GRATE_A, HYPER_GRATE_B*: the parameters for the beta prior of gain probability (default: 3,1).
  * *INIT_CONSERVE_RATE*: the initial conserved rate (default: 0.5).
  * *INIT_ACCE_RATE*: the initial accelerated rate (default: 1).
  * *CONSERVE_PRIOR_A*: the shape parameter for the gamma prior of conserved rate (default: 5).
  * *CONSERVE_PRIOR_B*: the scale parameter for the gamma prior of conserved rate (default: 0.04).
  * *ACCE_PRIOR_A*: the shape parameter for the gamma prior of accelerated rate (default: 10).
  * *ACCE_PRIOR_B*: the scale parameter for the gamma prior of accelerated rate (default: 0.2).
  * *RATE_OPT*: to avoid label switching of latent state, we provide options to restrict the range of accelerated rate and conserved rate. 0: no restriction on rates; 1: have lower and upper bound on accelerated and conserved rates respectively; 2: restrict accelerated rate to be larger than conserved rate. (default: 1)
  * *NLB*: lower bound for the accelerated rate. Only valid if *RATE_OPT* = 1. Default: 0.6.
  * *CUB*: upper bound for the conserved rate. Only valid if *RATE_OPT* = 1. Default: 1.

* **Control for MCMC and number of threads**: 
  * *BURNIN*: number of initial iterations to discard before equilibrium of the chain (default: 200). Should set to be larger.
  * *MCMC*: number of MCMC iterations (default: 800). Should set to be larger. 
  * *ADAPT_FREQ*: number of iterations to recompute acceptance rate of Metropolis-Hastings for adaptively adjusting the proposal variances for substitution rates (default: 500).
  * *SEED*: seed for random sampling (default: 5)
  * *SAMPLE_HYPER*: whether to sample hyperparameters. 0, fix hyperparameters; 1, sample (default: 0). Sampling hyperparameters is time-consuming, and not recommended. If sampling hyperparameters, the algorithm will only output the posterior of Z (latent state of each branch) under full model. 
  * *CHAIN*: Numer of iterations to sample hyper parameters. If not sampling hyperparameter, set it to 1 (default: 1).
  * *NUM_THREAD*: Number of threads to run the algorithm (default: 1).

