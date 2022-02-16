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

TODO:
- Split the Hu-etal-2019 folder into a separate repo within the PhyloAcc organization
- What is _config.yml? Seems to be related to building the website. Figure out where it needs to go
- Update docs
- Update links
- Remake the original PhyloAcc web page so that the link in the paper is still active. This page can just point to the new repo or page
- Make release for these changes and update conda recipe (meta.yaml) to point to it
- Fork bioconda and upload recipe
- Add a simple way to test the PhyloAcc-ST binary, maybe --version option and/or a very small test dataset
- Remove .py extension from interface (difficult for Windows)
- Resolve or transfer interface issues