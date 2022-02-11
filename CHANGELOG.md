Version 1.1.0, 02.11.2022
- Combining PhyloAcc and Python interface repos to facilitate codebase mergining in the future and to develop the conda package
- Moved all data files to Hu-etal-2019 to split into a separate repo later
- Updated Makefile to be used with conda build; the original Makefile can still be found in the src/PhyloACC-ST/ folder if we ever need to build from source that way
- Added meta.yaml, build.sh, and build.bat for conda building which will be split out later to our fork of bioconda
- Changed name of interface from 'phyloacc_interface.py' to 'phyloacc.py' and changed name of PhyloAcc binary to PhyloAcc-ST
- Changed default path for the PhyloAcc binary within the interface to point to PhyloAcc-ST to align with the name change above
- Changed licsene from MIT to GPL3
- Moved V2_GBGC to the src/PhyloAcc-ST/ folder





TODO:
- Split the Hu-etal-2019 folder into a separate repo within the PhyloAcc organization
- What is _config.yml? Seems to be related to building the website. Figure out where it needs to go
- Update docs
- Update links
- Remake the original PhyloAcc web page so that the link in the paper is still active. This page can just point to the new repo or page
- Make release for these changes and update conda recipe (meta.yaml) to point to it
- Fork bioconda and upload recipe
- Add a simple way to test the PhyloAcc-ST binary, maybe --version option and/or a very small test dataset
- Remove .py extension from interface
- Resolve or transfer interface issues