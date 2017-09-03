# PhyloAcc
This is a software to detect the shift pattern of DNA substitution rate of a genomic region, and identify genomic elements accelerated in some specific species from a set of conserved elements. The underlying model assumes a latent discrete state (Z) of relative substitution rate of each branch on the phylogeny which can be neutral, conserved and accelerated. For each genomic element, it will start from neutral rate at the common ancestor of the phylogeny, transit to conserved state and then reach accelerated state in some lineages. Our method utilizes adaptive collapsed Gibbs sampling to obtain the pattern of substitution rate shifts (posterior distribution of Z) as well as relative substitution rates of conserved and accelerated state. To identify DNA elements with accelerating on specific branches, we compared marginal likelihood under three models: null model (M0) where all branches are neutral or conserved; accelerated model (M1) in which branches leading to target species are accelerated; and full model (M2), with no constraints on latent states. Then we used two Bayes factors: between M1 and M0 (BF1) and between M1 and M2 (BF2) as criteria to identify DNA elements accelerated exclusively in target lineages.

## Getting Started
Some prelimenary inputs which might be generated from other softwares are required: (1) a Phylogeny in mod format. The file can be output from phyloFit in PHAST package, with the transition rate matrix of bases and branch lengths. Our model assumes that those are the substitution rates under neutral evoluation and will estimate conserved and accelerated rate relatived to the neutral rate. In our study, we used genome-wide four-fold degenerate sites to estimtate the rate matrix and branch lengths. (2) a multiple alignment file concatenating sequences of all input conserved elements in FASTA format, and (3) a bed file with the position of each individual element in the coordinate of concatenated alignment file (0-based).

We also need a parameter file, which contains the pathes for input files and output directory, information of species and parameters for MCMC. Please read [README2.md](https://github.com/xyz111131/PhyloAcc/blob/master/README2.md) for more detail. 

After running the algorithm, our method will output the posterior of latent state (Z) for each branch (indicating neutral, conserved or accelerated) for each element under each model in the files "*prefix*\_rate_postZ\_[0-2].txt" and the marginal loglikelihoods for each element are in the file "*prefix*_elem_lik.txt". The format of output files are explained in [README_OUTPUT.md](https://github.com/xyz111131/PhyloAcc/blob/master/README_OUTPUT.md).

## Prerequisites
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/):a numerical library for C and C++. 
* [Armadillo](http://arma.sourceforge.net/): C++ linear algebra library.
* [Open MPI](https://www.open-mpi.org/)(not required): for parallel computing.
## Installation

## Usage
PhyloAcc Simulation/param1-1.txt



