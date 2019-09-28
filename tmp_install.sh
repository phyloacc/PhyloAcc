#!/bin/bash

source activate phyloacc

export CONDA_ENV_PATH=$(conda info -e | grep "\*" | awk -F '[[:space:]]+' '{print $3}')
export LD_LIBRARY_PATH=${CONDA_ENV_PATH}/lib/:$LD_LIBRARY_PATH
export LD_RUN_PATH=${CONDA_ENV_PATH}/lib/:$LD_RUN_PATH

echo ${CONDA_ENV_PATH}

make
make install
