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

Starting with [v2.3.1](https://github.com/phyloacc/PhyloAcc/releases/tag/v2.3.1), there are now two ways to specify options for `phyloacc.py`:

1. Specifying them in a configuration file with `--config`
2. Specifying them in the command line

### Config file

Starting with [v2.3.1](https://github.com/phyloacc/PhyloAcc/releases/tag/v2.3.1), options for `phyloacc.py` can be specified either by the command line using the flags below or in a configuration file in [YAML format](https://en.wikipedia.org/wiki/YAML) and the using `--config <config file>` syntax.

We provide a configuration file template which lists all available options for `phyloacc.py` with no values specified: [Config file template](https://github.com/phyloacc/PhyloAcc/blob/main/phyloacc-cfg-template.yaml)

In this file, options are specified in key: value pairs, with the key being the option name on the left of the colon and the desired value being the to the right of the colon.

For example, I could run the program with the following command:

```bash
phyloacc.py -a alignment.fa -b loci.bed -m model.mod -t "species1;species2;species3"
```

**OR** I could specify these key value pairs in a config file `phyloacc-cfg.yaml`:

```yaml
aln_file: alignment.fa
bed_file: loci.bed
mod_file: model.mod
targets: species1;species2;species3
```

And then run:

```bash
phyloacc.py --config phyloacc-cfg.yaml
```

These are identical.

See the [Config file template](https://github.com/phyloacc/PhyloAcc/) for a full list of key names.

#### IMPORTANT: Options given via the command line take precedence over those given in the config file.

For instance, if I use the same config file above, but run the command:

```bash
phyloacc.py --config phyloacc-cfg.yaml -a other_alignment.fa
```

then `other_alignment.fa` will be used as the input alignment rather than `alignment.fa`.


If an option is not specified in either the command line or the config file, a default value will be used, which you can read about in the table below.

#### Boolean options

For boolean options in the config file (`theta_flag`, `dollo_flag`, `overwrite_flag`, `labeltree`, `summarize_flag`, `filter_alns`, `options_flag`, `depcheck`, `append_log_flag`, `info_flag`, `version_flag`, `quiet_flag`), type *true* or *false* (or leave blank for *false*). In other words:

```bash
phyloacc.py --config phyloacc-cfg.yaml --overwrite
```

is equivalent to setting

```yaml
overwrite_flag: true
```

in the config file and running

```bash
phyloacc.py --config phyloacc-cfg.yaml
```

In practice for many of these, using the command line option may be easier than changing the config file each time.

### Command line options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `--config [STRING]` | The path to a YAML formatted configuration file that specifies the program options (see above). Note that options specified on the command line take precedence of those specified in the config file. | Optional |
| `-a [FASTA FILE]` | An alignment file with all loci concatenated. `-b` must also be specified. Expected as FASTA format for now. | **One of `-a`/`-b` or `-d` is REQUIRED.** |
| `-b [BED FILE]` | A bed file with coordinates for the loci in the concatenated alignment file. `-a` must also be specified. | **One of `-a`/`-b` or `-d` is REQUIRED.**  |
| `-i [TEXT FILE]` | A text file with locus names, one per line, corresponding to regions in the input bed file. If provided, PhyloAcc will only be run on these loci. | Optional. **-a and -b must also be specified.**  |
| `-d [DIRECTORY]` | A directory containing individual alignment files for each locus. Expected as FASTA format for now. | **One of `-a`/`-b` or `-d` is REQUIRED.**  | 
| `-m [MOD FILE]` | A file with a background transition rate matrix and phylogenetic tree with branch lengths as output from phyloFit. | **REQUIRED.** |
| `-o [DIRECTORY]` | Desired output directory. This will be created for you if it doesn't exist. | phyloacc-[date]-[time]  |
| `-t "[STRING]"` | Tip labels in the input tree to be used as target species. Enter multiple labels separated by semi-colons (;). | **REQUIRED.** | 
| `-c "[STRING]"` | Tip labels in the input tree to be used as conserved species. Enter multiple labels separated by semi-colons (;). | Optional. Any species not specified in `-t` or `-g` will be inferred as conserved.  |
| `-g "[STRING]"` | Tip labels in the input tree to be used as outgroup species. Enter multiple labels separated by semi-colons (;). | Optional. |
| `-l [NEWICK FILE]` | A file containing a rooted, Newick formatted tree with the same topology as the species tree in the mod file (`-m`), but with branch lengths in coalescent units. | When the gene tree model is used, one of `-l` or `--theta` must be set. |
| `--theta` | Set this to add gene tree estimation with IQ-tree and species estimation with ASTRAL for estimation of the theta prior. Note that a species tree with branch lengths in units of substitutions per site is still required with `-m`. Also note that this may add substantial runtime to the pipeline. | When the gene tree model is used, one of `-l` or `--theta` must be set. |
| `-r [STRING]` | Determines which version of PhyloAcc will be used. gt: use the gene tree model for all loci, st: use the species tree model for all loci, adaptive: use the gene tree model on loci with many branches with low sCF and species tree model on all other loci. | st |
| `--dollo` | Set this to re-enable the Dollo assumption of non-reversibility from the original model. If set, once a branch is inferred to be in an accelerated state, it and its descendants will always be in an accelerated state. | Optional |
| `-n [INT]` | The number of processes that this script should use. | 1 |

### MCMC options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-mcmc [INT]` | The total number of steps in the Markov chain. | 1000 |
| `-burnin [INT]` | The number of steps to be discarded in the Markov chain as burnin. | 500 |
| `-chain [INT]` | The number of MCMC chains to run. | 1 |
| `-thin [INT]` | For the gene tree model, the number of MCMC steps between gene tree sampling. | 1 |

### sCF options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-scf [FLOAT]` | The value of sCF to consider as low for any given branch or locus. Must be between 0 and 1. | 0.5  |
| `-s [FLOAT]` | A value between 0 and 1. If provided, this proportion of branches must have sCF below `-scf` to be considered for the gene tree model. Otherwise, branch sCF values will be averaged for each locus. | Optional |

### Batching options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `-batch [INT]` | The number of loci to run per batch. | 50 |
| `-p [INT]` | The number of processes to use for each batch of PhyloAcc. | 1 |
| `-j [INT]` | The number of jobs (batches) to run in parallel. | 1 |

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
| `-st-path [STRING]` | The path to the PhyloAcc-ST binary. | PhyloAcc-ST |
| `-gt-path [STRING]` | The path to the PhyloAcc-ST binary. | PhyloAcc-GT |
| `-iqtree [STRING]` | When `--theta` is set, gene trees will be inferred with IQ-Tree. This specifies the path to the IQ-Tree executable. | iqtree |
| `-coal-cmd [STRING]` | When `--theta` is set, a species tree will be inferred with branch lengths coalescent branch units. This specifies the command to infer that tree from the gene trees estimated with IQ-Tree. Currently supported programs: ASTRAL | java -jar astral.jar |
| `-phyloacc "[STRING]"` | A catch-all option for other PhyloAcc parameters. Enter as a semi-colon delimited list of options: 'OPT1 value;OPT2 value' | Optional |
| `--options` | Print the full list of PhyloAcc options that can be specified with `-phyloacc` and exit. | Optional |

### Miscellaneous  options

| Option | Description | Default value |
| ------ | ----------- |---------------|
| `--labeltree` | Simply reads the tree from the input mod file (`-m`), labels the internal nodes, and exits.  | Optional |
| `--overwrite` | Set this to overwrite existing files in the specified output directory. | Optional |
| `--appendlog` | Set this to keep the old log file even if `--overwrite` is specified. New log information will instead be appended to the previous log file. | Optional |
| `--summarize` | Only generate the input summary plots and page. Do not write or overwrite batch job files.  | Optional |
| `--filter` | By default, alignments with 0 informative sites are removed from the analysis. Set this to also automatically filter alignments in which 50% of sites have many gaps (>50%), or 50% of sequences have many gaps (>50%) | Optional |
| `--info` | Print some meta information about the program and exit. No other options required. | Optional |
| `--depcheck` | Run this to check that all dependencies are installed at the provided path. No other options necessary. | Optional |
| `--version` | Simply print the version and exit. Can also be called as `-version`, `-v`, or `--v`. | Optional |
| `--quiet` | Set this flag to prevent PhyloAcc from reporting detailed information about each step. | Optional |
| `-h` | Print a help menu and exit. Can also be called as `--help`. | Optional |