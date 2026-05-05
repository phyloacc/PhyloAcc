TARGET_ST=PhyloAcc-ST
TARGET_GT=PhyloAcc-GT
# The name of the compiled binary

# make PREFIX=$CONDA_PREFIX
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/
# To make locally with the conda environment prefix

CXX=g++
# Which compiler to use
# Note: g++ 5.4 resulted in several errors while compiling: SRC/bpp_c2.cpp:345:12: error: ‘::isnan’ has not been declared
# Require g++ 7+

$(info $$CXX is [${CXX}])
# Report the compiler used

$(info $$PREFIX is [${PREFIX}])
# Report the PATH prefix

GSL_INCLUDE=${PREFIX}/include/
GSL_LIB=${PREFIX}/lib/
$(info $$GSL_INCLUDE is [${GSL_INCLUDE}])
$(info $$GSL_LIB is [${GSL_LIB}])
# GSL paths with the conda environment prefix

CFLAGS=-Wall -g -O2 -std=c++14
LDFLAGS=-lgsl -lm -lgslcblas -larmadillo -fopenmp -Wl,-rpath,$(GSL_LIB)
# Options for the g++ commands

############

SRC_DIR_COMMON=src/PhyloAcc-common/
SRCS_COMMON=$(SRC_DIR_COMMON)/*.cpp
INCLUDES_COMMON=$(SRC_DIR_COMMON)/*.h
COMMON_INCLUDE=-I$(SRC_DIR_COMMON)

SRC_DIR_ST=src/$(TARGET_ST)/
# Shared implementations now live in src/PhyloAcc-common/.
# The legacy ST-local profile/newick/utils copies are intentionally excluded
# from the build and should not be treated as the active implementation path.
SRCS_ST=$(wildcard $(SRC_DIR_ST)*.cpp) $(SRCS_COMMON)
INCLUDES_ST=$(wildcard $(SRC_DIR_ST)*.h) $(wildcard $(SRC_DIR_ST)*.hpp) $(INCLUDES_COMMON)
ST_INCLUDE=-I$(SRC_DIR_ST)
# Locations of files to compile

$(TARGET_ST): $(SRCS_ST) $(INCLUDES_ST)
	$(CXX) $(CFLAGS) $(COMMON_INCLUDE) $(ST_INCLUDE) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(SRCS_ST) -o $(TARGET_ST) $(LDFLAGS)
# g++ commands for each file
# Species tree version
############

SRC_DIR_GT=src/$(TARGET_GT)/
# Shared implementations now live in src/PhyloAcc-common/.
# The legacy GT-local profile/newick/utils copies are intentionally excluded
# from the build and should not be treated as the active implementation path.
SRCS_GT=$(wildcard $(SRC_DIR_GT)*.cpp) $(SRCS_COMMON)
INCLUDES_GT=$(wildcard $(SRC_DIR_GT)*.h) $(wildcard $(SRC_DIR_GT)*.hpp) $(INCLUDES_COMMON)
GT_INCLUDE=-I$(SRC_DIR_GT)
# Locations of files to compile

$(TARGET_GT): $(SRCS_GT) $(INCLUDES_GT)
	$(CXX) $(CFLAGS) $(COMMON_INCLUDE) $(GT_INCLUDE) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(SRCS_GT) -o $(TARGET_GT) $(LDFLAGS)
# g++ commands for each file
# Gene tree version
############

.PHONY: install
install: $(TARGET_ST) $(TARGET_GT)
	cp $(TARGET_ST) $(PREFIX)/bin/$(TARGET_ST)
	cp $(TARGET_GT) $(PREFIX)/bin/$(TARGET_GT)
# Command to install by moving binary

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/bin/$(TARGET_ST)
	rm -f $(PREFIX)/bin/$(TARGET_GT)
# Command to uninstall by removing binary

.PHONY: clean
clean:
	rm -f *.o *~ $(TARGET_ST)
	rm -f *.o *~ $(TARGET_GT)
# Command to remove all compiled files to make a clean install

.PHONY: cpp-tests
cpp-tests:
	$(MAKE) -C tests/cpp
