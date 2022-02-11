TARGET=PhyloAcc-ST
# The name of the compiled binary

ifeq ($(OS),Windows_NT)
	CXX=g++
	GSL_INCLUDE=${LIBRARY_INC}
	GSL_LIB=${LIBRARY_LIB}
# GSL paths with the conda environment prefix
else
	CXX=g++-7
	GSL_INCLUDE=${PREFIX}/include/
	GSL_LIB=${PREFIX}/lib/
# GSL paths with the conda environment prefix
endif

# ifeq ($(shell uname),Linux Darwin)
# 	CXX=g++-7
# else
# 	CXX=g++
# endif
#CXX=g++-7
# Which compiler to use.
# Note: g++ 5.4 resulted in several errors while compiling: SRC/bpp_c2.cpp:345:12: error: ‘::isnan’ has not been declared
# Switched to g++-7

$(info $$CXX is [${CXX}])
# Report the compiler used

$(info $$PREFIX is [${PREFIX}])
$(info $$GSL_INCLUDE is [${GSL_INCLUDE}])
$(info $$GSL_LIB is [${GSL_LIB}])

CFLAGS=-Wall -g -O2 -std=c++11
LDFLAGS=-lgsl -lm -lgslcblas -larmadillo -fopenmp
# Options for the g++ commands

SRC_DIR=src/$(TARGET)/
SRCS=$(SRC_DIR)/*.cpp
INCLUDES=$(SRC_DIR)/*.h $(SRC_DIR)/*.hpp
# Locations of files to compile

#export LD_LIBRARY_PATH=${PREFIX}/lib/
#export LD_RUN_PATH=${PREFIX}/lib/
# LD paths for the dependencies installed through conda (armadillo, gsl, blas, etc.)
# e.g.:
# libcblas.so.3, needed by /home/gregg/anaconda3/envs/phyloacc/lib//libgsl.so, not found
# However, these don't appear to be needed with conda

$(TARGET): $(SRCS) $(INCLUDES)
	$(CXX) $(CFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(SRCS) -o $(TARGET) $(LDFLAGS)
# g++ commands for each file

.PHONY: install
install: $(TARGET)
ifeq ($(OS),Windows_NT)
	cp $< $(LIBRARY_BIN)/$(TARGET)
# GSL paths with the conda environment prefix
else
	cp $< $(PREFIX)/bin/$(TARGET)
endif
# Command to install by moving binary

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)
# Command to uninstall by removing binary

.PHONY: clean
clean:
	rm -f *.o *~ $(TARGET)
# Command to remove all compiled files to make a clean install



# Install notes:
# conda install -c conda-forge armadillo -> includes: conda install lapack
# conda install -c conda-forge gsl
# no package called: conda install atlas
#
# also had to: sudo apt install make
# sudo apt update
# sudo apt install g++ gdb make ninja-build rsync zip






