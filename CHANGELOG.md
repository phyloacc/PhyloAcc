Version 1.1.0, 02.11.2022
- Combining PhyloAcc and Python interface repos to facilitate codebase mergining in the future and to develop the conda package
- Moved all data files to Hu-etal-2019 to split into a separate repo later
- Updated Makefile to be used with conda build; the original Makefile can still be found in the src/PhyloACC-ST/ folder if we ever need to build from source that way
- Added meta.yaml, build.sh, and build.bat for conda building which will be split out later to our fork of bioconda
- Changed name of interface from 'phyloacc_interface.py' to 'phyloacc.py' and changed name of PhyloAcc binary to PhyloAcc-ST
- Changed default path for the PhyloAcc binary within the interface to point to PhyloAcc-ST to align with the name change above
- Changed licsene from MIT to GPL3
- Moved V2_GBGC to the src/PhyloAcc-ST/ folder

02.15.2022, interface
- Fixed coloring of branches to include internal branches in the summary page
- Added group input specification by internal nodes
- Added option to only re-run the summary and plot generation and ignore/don't write or overwrite the job files
- Fixed issue where log info wasn't actually written to the log
- Added an option to append to a previous log rather than overwrite

02.16.2022, interface
- Added a check for duplicate labels in input groups
- Removed --plot option and changed --plotonly option to --summarize; plots are now always generated and --summarize indicates job files should not be overwritten or generated
- Fixed error codes in opt_parse

02.16.2022
- Moved the V2_GBGC source from within the PhyloAcc-ST dir to its own dir in src called PhyloAcc-ST-GBGC
- Removed the simulation data from the GBGC folder and put it in the Hu-etal-2019 folder under the GBGC subfolder
- Added test folder with original test data from the README (500 simulated elements on the ratite tree, but with the ID file only 10 elements are run)

02.18.2022
- Split the Hu-etal-2019 subfolder into its own repository
- Move the test data to a new repository, PhyloAcc-test-data
