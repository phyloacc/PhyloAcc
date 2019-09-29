#!/bin/bash

# Install gcc, conda, and git.
# We suggest you use brew to install cmake and gcc
# brew install gcc
# brew install git

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
# Note: in this way, we install them only for this env.
conda install gsl
conda install -c conda-forge openmp

# get conda env path
# Note: please check this value for yourself.
export CONDA_ENV_PATH=$(conda info -e | grep "\*" | awk -F '[[:space:]]+' '{print $3}')
# LD_LIBRARY_PARH is used for the program to find the libraries when it runs.
export LD_LIBRARY_PATH=${CONDA_ENV_PATH}/lib/:$LD_LIBRARY_PATH
# The two variables are used to tell make how to find the gsl, openmp when compiling the program.
# Note we install armadillo with brew. In fact, we test conda install armadillo, but it has some problems
# when compile. Help is needed.
export CONDA_ENV_INCLUDE=${CONDA_ENV_PATH}/include/
export CONDA_ENV_LIB=${CONDA_ENV_PATH}/lib/
# export LD_RUN_PATH=${CONDA_ENV_PATH}/lib/:$LD_RUN_PATH

# Even in macOS, we use g++ to compile (not clang).
# In makefile, we use g++, so if you have no g++ in your local bin path,
# you can use soft link below. (Ideally, after brewinng install gcc, it will
# link gcc to the installed dir by brew, i.e., /usr/local/Cellar/gcc).
ln -s /usr/local/Cellar/gcc/9.2.0/bin/g++-9 /usr/local/bin/g++

# Now, compile the program.
make
# Done. A local executable file called PhyloAcc should in the current directory.
