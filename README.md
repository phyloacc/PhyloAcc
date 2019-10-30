# PhyloAcc
This is a software to detect the shift pattern of DNA substitution rate of a genomic region and identify genomic elements accelerated in some specific species from a set of conserved elements. The underlying model assumes a latent discrete state (Z) of relative substitution rate of each branch on the phylogeny which can be neutral, conserved and accelerated. For each genomic element, it will start from a neutral or conserved state at the common ancestor of the phylogeny, transit to conserved state if not yet being conserved and then reach an accelerated state in some lineages. Our method utilizes adaptive collapsed Gibbs sampling to obtain the pattern of substitution rate shifts (posterior distribution of Z) as well as relative substitution rates of conserved and accelerated state. To identify DNA elements with accelerating on specific branches, we compare marginal likelihoods under three models: null model (M0) where all branches are neutral or conserved; accelerated model (M1) in which branches leading to target species are accelerated; and full model (M2), with no constraints on latent states. Then we use two Bayes factors: between M1 and M0 (BF1) and between M1 and M2 (BF2) as criteria to identify DNA elements accelerated exclusively in target lineages.

## Getting Started
Some preliminary inputs which might be generated from other software are required: (1) a Phylogeny in mod format. The file can be output from phyloFit in PHAST package, with the transition rate matrix of bases and branch lengths. Our model assumes that these branch lengths represent the expected number of substitutions in the background state and will estimate the conserved and accelerated rate relative to the background rate. In our study, we used genome-wide four-fold degenerate sites to estimate the rate matrix and branch lengths. (2) a multiple alignment files concatenating sequences of all input conserved elements in FASTA format, and (3) a bed file with the position of each individual element in the coordinate of concatenated alignment file (0-based).


We also need a parameter file, which contains the paths for input files and output directory, information of species and parameters for MCMC. Please read [README_PARAMETER.md](https://github.com/xyz111131/PhyloAcc/blob/master/README_PARAMETER.md) for more detail. 


After running the algorithm, our method will output the posterior of latent state (Z) for each branch (indicating background, conserved or accelerated) for each element under each model in the files "*prefix*\_rate_postZ\_M[0-2].txt" and the marginal loglikelihoods for each element are in the file "*prefix*_elem_lik.txt". The format of output files is explained in [README_OUTPUT.md](https://github.com/xyz111131/PhyloAcc/blob/master/README_OUTPUT.md).

## Prerequisites
* [GCC](https://gcc.gnu.org/): You might need latest GCC (version 7) supporting openmp. If you are using Mac, you could use brew to (re)install gcc. 
```bash
brew update ## update the formulae and Homebrew itself, if your brew is out-dated
brew install gcc
```
* [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/): a numerical library for C and C++. PhyloAcc has been
  tested with version 2.4 of GSL.
* [Armadillo](http://arma.sourceforge.net/): C++ linear algebra library. You could install Armadillo following the steps on its website. For Linux, before installing Armadillo, you need to install CMAKE, LAPACK, BLAS (OPENBLAS) and ATLAS, along with the corresponding development/header files. PhyloAcc has been tested with version 7.800.2.
For Mac, you could use brew (tested and Recommended): 
```bash
brew install homebrew/science/armadillo
```
* [Open MP](http://www.openmp.org/): for parallel computing. 
* To use the R functions to plot,  please install [Rstudio](https://www.rstudio.com/) with the current version of R (>=3.3.2) and install seqinr, ggplot2, reshape2, ape packages.  

To install these softwares on Mac, you could follow the procedures below:
1) make a new conda environment called "PhyloAcc" for example, and activate the environment
2) conda install gsl
3) conda install atlas
4) conda install lapack
5) conda install armadillo
6) edit your .bashrc in your home directory:
   ```bash
    export LD_LIBRARY_PATH=/user/miniconda3/envs/PhyloAcc/lib/:$LD_LIBRARY_PATH
    export LD_RUN_PATH=/user/miniconda3/envs/PhyloAcc/lib/:$LD_RUN_PATH
    ```
    or other path of the conda environment
7) source .bashrc
8) change to g++-9 in the PhyloAcc Makefile if you are using GCC version 9
9) edit the location of GSL lib and include in the Makefile: e.g. /user/miniconda3/envs/PhyloAcc
10) edit *PREFIX* in the Makefile for the installation directory of the conda env

*(credited to Wei Gordon)*

## Build on Linux or Mac
Run:
```bash
make
```
in PhyloAcc directory to generate the 'PhyloAcc' executable.

## Installation
Run:
```bash
sudo make install
```
to install in default path /usr/local/bin, and 
```bash
sudo make uninstall
```
to uninstall.

For our extended version modeling GC-based gene conversion (gBGC) effect, please go to [V2_GBGC/](https://github.com/xyz111131/PhyloAcc/blob/master/V2_GBGC) and make & install under that directory.

## Usage
Try this in PhyloAcc directory as a test, which will run simulated elements on ratite phylogeny:
```bash
mkdir Simulation_ratite/result_tmp
./PhyloAcc Simulation_ratite/param2-1-test.txt
```
or this after installation:
```bash
PhyloAcc Simulation_ratite/param2-1-test.txt
```
For testing propose, it will only run the first 10 elements of simulated data from Simulation_ratite/simu_500_200_diffr_2-1.* and output to Simulation_ratite/result_tmp/. To run all the elements and get results in Simulation_ratite/result_phyloAcc/, you could run:
```bash
./PhyloAcc Simulation_ratite/param2-6.txt
```
To run your own data, please change the paths in your parameter file.

To run the model including gBGC,
```bash
cd V2_GBGC
mkdir Simulation/result_tmp
./PhyloAcc_gBGC paramGC-0.txt
```
under *V2_GBGC/*. It will output to *V2_GBGC/Simulation/result_tmp/* (which will be the same as in V2_GBGC/Simulation/result/).

There are several R scripts available in [R/](https://github.com/xyz111131/PhyloAcc/blob/master/R) which read the output from PhyloAcc and generate plots in the main paper (e.g. "scaled" phylogenetic tree and sequence alignment for one element). Please read [plot.html](https://xyz111131.github.io/PhyloAcc/R/plot.html) and run [plot.Rmd](https://github.com/xyz111131/PhyloAcc/blob/master/R/plot.Rmd) for detail. R scripts to generate simulated DNA sequences are also in [R/](https://github.com/xyz111131/PhyloAcc/blob/master/R), please see [Simulation.md](https://github.com/xyz111131/PhyloAcc/blob/master/Simulation.md) for more detail. 

## Data
The species names and phylogenetic trees used in the PhyloAcc manuscript are in [Data/](https://github.com/xyz111131/PhyloAcc/blob/master/Data/). The simulated sequences and results are in [Simulation_mammal/](https://github.com/xyz111131/PhyloAcc/blob/master/Simulation_mammal/) and [Simulation_ratite/](https://github.com/xyz111131/PhyloAcc/blob/master/Simulation_ratite/). The results for ratite and mammalian CNEEs in [mammal_result/](https://github.com/xyz111131/PhyloAcc/blob/master/mammal_result/) and [ratite_result/](https://github.com/xyz111131/PhyloAcc/blob/master/ratite_result/). 

## Trouble Shooting
After running phyloFit for background sequences to get the rate matrix and branch lengths, you might need to name all ancestral nodes in the tree. You could use tree_doctor in phast (http://compgen.cshl.edu/phast/help-pages/tree_doctor.txt): 
``` bash
tree_doctor --name-ancestors input.mod > output.mod
```

## Reference
Zhirui Hu, Timothy B Sackton, Scott V. Edwards, Jun S. Liu: **A hierarchical Bayesian model for detecting convergent rate changes of conserved noncoding elements on phylogenetic trees**, *bioRxiv*, 2018.
https://www.biorxiv.org/content/early/2018/02/22/260745.1
