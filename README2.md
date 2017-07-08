# Parameter file
In the parameter file, each parameter is specified in a line with the parameter name at the beginning followed by paremeter value. The parameters are: 
* Input and output: 
  * **PHYTREE_FILE**: the path of Phylogeny (.mod)  
  * **SEG_FILE**: the path of bed file for genomic regions
  * **ALIGN_FILE**: the path of multiple alignment file (.fasta)
  * **RESULT_PREFIX**: the prefix for output files
  * **BURNIN**: number of initial iterations to discard before equilibrium of the chain (default: 500)

* Control for MCMC: 
MCMC: number of MCMC iterations (default: )
CHAIN: 1
TARGETSPECIES: species of interest. E.g. species potentially lost conservation or with convergent phenotype changes.
OUTGROUP: outgroup species of the phylogeny. These species are not considered to be accelerated in our model. 
CONSERVE: species assumed to be mostly conserved. The algorithm will filter out elements with alignment gaps in more than 50% of the conserved species. 
REF: galGal
INDEL 0
GAPCHAR -
