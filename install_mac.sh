#!/bin/bash

# Install cmake, gcc, conda, and git.
# We suggest you use brew to install cmake and gcc
brew install cmake
brew install gcc

# You can install install conda from their website.
# We use conda as the basic vitural enviorment for the package.

# Install third party package: armadillo, openmp
# Note: when install armadillo, the system might say warning,
# which should be OK.
brew install armadillo

# Git clone the repo to your local path.
# The command below will download the repo in your current path.
git clone https://github.com/xyz111131/PhyloAcc.git
# Enter this directory
cd PhyloAcc

# Create the conda env for phyloacc
conda create -n phyloacc
# enter that env
source activate phyloacc

# Install gsl, openmp
# Note: in this way, we install gsl only for this env.
conda install gsl
conda install -c conda-forge openmp

# get conda env path
export CONDA_ENV_PATH=$(conda info -e | grep "\*" | awk -F '[[:space:]]+' '{print $3}')
export LD_LIBRARY_PATH=${CONDA_ENV_PATH}/lib/:$LD_LIBRARY_PATH
export CONDA_ENV_INCLUDE=${CONDA_ENV_PATH}/include/
export CONDA_ENV_LIB=${CONDA_ENV_PATH}/lib/
# export LD_RUN_PATH=${CONDA_ENV_PATH}/lib/:$LD_RUN_PATH

# Check if the variable below has value by print it.
echo ${CONDA_ENV_PATH}

# Even in macOS, we use g++ to compile (not clnag).
# In makefile, we use g++, so if you have no g++ in your local bin path,
# you can use soft link below. (Ideally, after brewinng install gcc, it will
# link gcc to the brew installed dir, i.e., /usr/local/Cellar/gcc).
ln -s /usr/local/Cellar/gcc/9.2.0/bin/g++-9 /usr/local/bin/g++

# Now, compile the program.
make
# Done. A local executable file called PhyloAcc should in the current directory.
