Version 2.4.3, 03.25.2025
- Fixed how tree labels are read in GT
- Switched from bash script to Python functions to label ASTRAL tree internal nodes, which required adding the various hidden `-lt` options to `phyloacc.py`

Version 2.4.0, 09.26.2024
- Added `--local` option to generate a snakemake command that does not use a cluster, intended for testing purposes only.
- Fixed the `--version`, `--quiet`, and `--options` flags to work with the config file.
- Fixed spacing when string is longer than the pad.
- Added `--testcmd` to also print a direct PhyloAcc command at the end of the interface for testing purposes.

Version 2.3.4, 09.24.2024
- Reverted the element ID indexing in the interface from 1 to 0 to match the C++ code
- Changed how the interface uses regex to read trees so it uses `r`aw strings
- Fixed bug in which gene trees were still being inferred with the ST model

Version 2.3.0
- Check for infs in `phyloacc_post.py`
- Added error checking for unlabeled trees and alignments with labels that don't match the tree
- Added capability to specify input options in a config file with `--config`
- Provided template config file (`phyloacc-cfg.yml`)
- Added the `--filter` option to filter out alignments with too many missing sites
- Added hidden `--debug-aln` option to stop the program after reading the alignments

Version 2.2.0, 04.13.2023
- Added `--nophyloacc` option that prevents execution of the PhyloAcc rules in the snakemake workflow, useful for debugging or just running `--theta`
- Internally, switched the number of informative sites required for a locus to be used in `--theta` estimation to be a param, maybe user option later
- Added `--dollo` option which sets the PhyloAcc `HYPER_LRATE2_A` parameter to 0 to use the Dollo assumption from the original model
- Implemented (better) handling of Warnings
- Added hidden `--dev` option to automatically set paths when testing things locally
- Added hidden `-inf-frac-theta` option to control the fraction of informative sites needed for a locus to be used in `--theta` estimation (for development)
- Fixed bug in which options that require floats as input (e.g. `-scf`) wouldn't be set if the input value was 0.0

04.06.2023
- Fixed ids dir when `--theta` is set
- Added `-iqtree-path` and `-coal-cmd` options to specify programs to estimate branch lengths in coalescent units when `--theta` is set
- Fixed bug where `--labeltree` added NA labels to nodes in addition to <#> nodes when input tree is unlabeled
- Added error message if tree in .mod file does not have labels on the internal nodes telling the user that PhyloAcc requires a labeled tree as input
- `--labeltree` now also prints a tree labeled in the format of [descendant 1]-[descendant 2], where descendants 1 and 2 are the first of the alphabetically sorted tips that descend from that node

Version 1.1.0, 02.11.2022
- Combining PhyloAcc and Python interface repos to facilitate codebase mergining in the future and to develop the conda package
- Moved all data files to `Hu-etal-2019` to split into a separate repo later
- Updated Makefile to be used with conda build; the original Makefile can still be found in the `src/PhyloACC-ST/` folder if we ever need to build from source that way
- Added `meta.yaml`, `build.sh`, and `build.bat` for conda building which will be split out later to our fork of bioconda
- Changed name of interface from `phyloacc_interface.py` to `phyloacc.py` and changed name of `PhyloAcc` binary to `PhyloAcc-ST`
- Changed default path for the PhyloAcc binary within the interface to point to `PhyloAcc-ST` to align with the name change above
- Changed licsene from MIT to GPL3
- Moved `V2_GBGC` to the `src/PhyloAcc-ST/` folder

02.15.2022, interface
- Fixed coloring of branches to include internal branches in the summary page
- Added group input specification by internal nodes
- Added option to only re-run the summary and plot generation and ignore/don't write or overwrite the job files
- Fixed issue where log info wasn't actually written to the log
- Added an option to append to a previous log rather than overwrite

02.16.2022, interface
- Added a check for duplicate labels in input groups
- Removed `--plot` option and changed `--plotonly` option to `--summarize`; plots are now always generated and `--summarize` indicates job files should not be overwritten or generated
- Fixed error codes in opt_parse

02.16.2022
- Moved the `V2_GBGC` source from within the `PhyloAcc-ST` dir to its own dir in src called `PhyloAcc-ST-GBGC`
- Removed the simulation data from the `GBGC` folder and put it in the `Hu-etal-2019` folder under the `GBGC` subfolder
- Added `test` folder with original test data from the README (500 simulated elements on the ratite tree, but with the ID file only 10 elements are run)

02.18.2022
- Split the `Hu-etal-2019` subfolder into its own repository
- Move the test data to a new repository, `PhyloAcc-test-data`

03.07.2022
- Remove `rcParams["xtick.labelcolor"]` and `rcParams["ytick.labelcolor"]` since they don't seem to work for `matplotlib-base`
- Remove `--plot` option from post script; plots are now generated by default
- Change version string

07.15.2022
- Implemented tree parsing as a class rather than functions
- Refactored sCF code
- Organized option parsing and main code
- Added `-scf` and `-s` options for sCF user cutoffs
- Updated README and moved version 1 READMEs to `docs/`

10.21.2022
- More testing of tree class
- Added `thin` as a user option
- Split the HTML functions out of `plot.py` into `html.py`
- Added new plots and output tables for `phyloacc_post.py`
- Fixed default path of config file in PhyloAcc-GT `main.cpp`, now exits rather than trying a non-existant file
