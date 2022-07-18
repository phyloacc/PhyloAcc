# <p align="center"><a href="https://phyloacc.github.io/" target="_blank"><img align="center" width="360" height="210" src="https://phyloacc.github.io/img/logo2.png"></a></p>

## PhyloAcc: Substitution rate estimate and rate shift inference along branches of a phylogeny

[![Install](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://phyloacc.github.io/install.html)
[![OS](https://anaconda.org/bioconda/phyloacc/badges/platforms.svg)](https://phyloacc.github.io/index.html)
[![Version](https://img.shields.io/conda/vn/bioconda/phyloacc?label=version)](https://bioconda.github.io/recipes/phyloacc/README.html)
[![Release Date](https://anaconda.org/bioconda/phyloacc/badges/latest_release_date.svg)](https://bioconda.github.io/recipes/phyloacc/README.html)
[![Downloads](https://img.shields.io/conda/dn/bioconda/phyloacc.svg?style=flat)](https://bioconda.github.io/recipes/phyloacc/README.html)
[![Commits](https://img.shields.io/github/commits-since/phyloacc/PhyloAcc/v2.0.0)](https://github.com/phyloacc/PhyloAcc/commits/main)
[![License](https://anaconda.org/bioconda/phyloacc/badges/license.svg)](https://github.com/phyloacc/PhyloAcc/blob/main/LICENSE.md)

# Authors
### Han Yan, Zhirui Hu, Gregg Thomas, Scott Edwards, Jun Liu, Tim Sackton

---
</br>

# Please see [the PhyloAcc website](https://phyloacc.github.io/) for full installation and usage instructions. What follows is only basic info.
</br>

# Basic info

## Installation

If conda/bioconda is already set up on your system, you can install PhyloAcc with a single command:

```bash
conda install phyloacc
```

For more detailed instructions and troubleshooting, see [the Installation page](https://phyloacc.github.io/install.html). If you have other questions or trouble let us know with [an issue](https://github.com/phyloacc/PhyloAcc/issues).

## Usage

For more detailed information and example commands, see [the README on the PhyloAcc website](https://phyloacc.github.io/readme.html)

### Input

* A species tree and background substitution rate estimates in `.mod` file format from [phyloFit](http://compgen.cshl.edu/phast/phyloFit-tutorial.php). [Example file](https://github.com/phyloacc/PhyloAcc-test-data/blob/main/ratite.mod)

* Either a single concatenated alignment in FASTA format with all elements to test ([Example file](https://github.com/phyloacc/PhyloAcc-test-data/blob/main/simu_500_200_diffr_2-1.fa)) AND a bed file that designates the coordinates of the elements ([Example file](https://github.com/phyloacc/PhyloAcc-test-data/blob/main/simu_500_200_diffr_2-1.bed); note that only the first 3 columns are required) __OR__ a directory containing separate alignments for each element in FASTA format.

* For the gene tree model, a species tree with the same topology as the one in the `.mod` file, but with branch lengths in coalescent units. If this isn't easily available one can be estimated from the input elements with the `--theta` option, though this will increase runtime and using these elements may not result in the most accurate branch length estimates.

## Options

**Use `phyloacc.py -h` to print out a help menu listing all the options.**

### Input/output options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-a [FASTA FILE]` | An alignment file with all loci concatenated. `-b` must also be specified. Expected as FASTA format for now. | **One of `-a`/`-b` or `-d` is REQUIRED.** |
| `-b [BED FILE]` | A bed file with coordinates for the loci in the concatenated alignment file. `-a` must also be specified. | **One of `-a`/`-b` or `-d` is REQUIRED.**  |
| `-i [TEXT FILE]` | A text file with locus names, one per line, corresponding to regions in the input bed file. If provided, PhyloAcc will only be run on these loci. | Optional. **-a and -b must also be specified.**  |
| `-d [DIRECTORY]` | A directory containing individual alignment files for each locus. Expected as FASTA format for now. | **One of `-a`/`-b` or `-d` is REQUIRED.**  | 
| `-m [MOD FILE]` | A file with a background transition rate matrix and phylogenetic tree with branch lengths as output from phyloFit. |**REQUIRED.** |
| `-o [DIRECTORY]` | Desired output directory. This will be created for you if it doesn't exist. | phyloacc-[date]-[time]  |
| `-t "[STRING]"` | Tip labels in the input tree to be used as target species. Enter multiple labels separated by semi-colons (;). |**REQUIRED.** | 
| `-c "[STRING]"` | Tip labels in the input tree to be used as conserved species. Enter multiple labels separated by semi-colons (;). | Optional. Any species not specified in `-t` or `-g` will be inferred as conserved.  |
| `-g "[STRING]"` | Tip labels in the input tree to be used as outgroup species. Enter multiple labels separated by semi-colons (;). | Optional. |
| `-l [NEWICK FILE]` | A file containing a rooted, Newick formatted tree with the same topology as the species tree in the mod file (`-m`), but with branch lengths in coalescent units. | When the gene tree model is used, one of `-l` or `--theta` must be set. |
| `--theta` | Set this to add gene tree estimation with IQ-tree and species estimation with ASTRAL for estimation of the theta prior. Note that a species tree with branch lengths in units of substitutions per site is still required with `-m`. Also note that this may add substantial runtime to the pipeline. | When the gene tree model is used, one of `-l` or `--theta` must be set. |
| `-r [STRING]` | Determines which version of PhyloAcc will be used. gt: use the gene tree model for all loci, st: use the species tree model for all loci, adaptive: use the gene tree model on loci with many branches with low sCF and species tree model on all other loci. | st  |
| `-n [INT]` |  	The number of processes that this script should use. | 1 |

### MCMC options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-burnin [INT]` | The number of steps to be discarded in the Markov chain as burnin. | 500 |
| `-mcmc [INT]` | The total number of steps in the Markov chain. | 1000 |

### sCF options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-scf [FLOAT]` | The value of sCF to consider as low for any given branch or locus. Must be between 0 and 1. | 0.5  |
| `-s [FLOAT]` | A value between 0 and 1. If provided, this proportion of branches must have sCF below `-scf` to be considered for the gene tree model. Otherwise, branch sCF values will be averaged for each locus. | Optional |

### Batching options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-p [INT]` | The number of processes to use for each batch of PhyloAcc. | 1 |
| `-j [INT]` | The number of jobs (batches) to run in parallel. | 1 |
| `-batch [INT]` | The number of loci to run per batch. | 50 |

### Cluster options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-part "[STRING]"` | The SLURM partition or list of partitions (separated by commas) on which to run PhyloAcc jobs. | **REQUIRED** |
| `-nodes [INT]` | The number of nodes on the specified partition to submit jobs to. | 1 |
| `-mem [INT]` | The max memory for each job in GB. | 4 |
| `-time [INT]` | The time in hours to give each job. | 1 |

### Other PhyloAcc options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-path [STRING]` | The path to the PhyloAcc binary. | PhyloAcc |
| `-phyloacc "[STRING]"` | A catch-all option for other PhyloAcc parameters. Enter as a semi-colon delimited list of options: 'OPT1 value;OPT2 value' | Optional |
| `--options` | Print the full list of PhyloAcc options that can be specified with `-phyloacc` and exit. | Optional |

### Miscellaneous  options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `--labeltree` | Simply reads the tree from the input mod file (`-m`), labels the internal nodes, and exits.  | Optional |
| `--overwrite` | Set this to overwrite existing files in the specified output directory. | Optional |
| `--appendlog` | Set this to keep the old log file even if `--overwrite` is specified. New log information will instead be appended to the previous log file. | Optional |
| `--summarize` | Only generate the input summary plots and page. Do not write or overwrite batch job files.  | Optional |
| `--info` | Print some meta information about the program and exit. No other options required. | Optional |
| `--depcheck` | Run this to check that all dependencies are installed at the provided path. No other options necessary. | Optional |
| `--version` | Simply print the version and exit. Can also be called as `-version`, `-v`, or `--v`. | Optional |
| `--qiuet` | Set this flag to prevent PhyloAcc from reporting detailed information about each step. | Optional |
| `-h` | Print a help menu and exit. Can also be called as `--help`. | Optional |