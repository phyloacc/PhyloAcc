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
LDFLAGS=-lgsl -lm -lgslcblas -larmadillo -fopenmp
# Options for the g++ commands

############

SRC_DIR_ST=src/$(TARGET_ST)/
SRCS_ST=$(SRC_DIR_ST)/*.cpp
INCLUDES_ST=$(SRC_DIR_ST)/*.h $(SRC_DIR_ST)/*.hpp
# Locations of files to compile

$(TARGET_ST): $(SRCS_ST) $(INCLUDES_ST)
	$(CXX) $(CFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(SRCS_ST) -o $(TARGET_ST) $(LDFLAGS)
# g++ commands for each file
# Species tree version
############

SRC_DIR_GT=src/$(TARGET_GT)/
SRCS_GT=$(SRC_DIR_GT)/*.cpp
INCLUDES_GT=$(SRC_DIR_GT)/*.h $(SRC_DIR_GT)/*.hpp
# Locations of files to compile

$(TARGET_GT): $(SRCS_GT) $(INCLUDES_GT)
	$(CXX) $(CFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(SRCS_GT) -o $(TARGET_GT) $(LDFLAGS)
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
