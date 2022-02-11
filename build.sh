#!/bin/bash

make 
make install
cp src/interface/phyloacc_interface.py ${PREFIX}/bin/.
mkdir -p ${SP_DIR}
cp -R src/interface/phyloacc_lib ${SP_DIR}/.
