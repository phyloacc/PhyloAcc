TARGET=PhyloAcc-ST
# The name of the compiled binary

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

CFLAGS=-Wall -g -O2 -std=c++11
LDFLAGS=-lgsl -lm -lgslcblas -larmadillo -fopenmp
# Options for the g++ commands

SRC_DIR=src/$(TARGET)/
SRCS=$(SRC_DIR)/*.cpp
INCLUDES=$(SRC_DIR)/*.h $(SRC_DIR)/*.hpp
# Locations of files to compile

$(TARGET): $(SRCS) $(INCLUDES)
	$(CXX) $(CFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(SRCS) -o $(TARGET) $(LDFLAGS)
# g++ commands for each file

.PHONY: install
install: $(TARGET)
	cp $< $(PREFIX)/bin/$(TARGET)
# Command to install by moving binary

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)
# Command to uninstall by removing binary

.PHONY: clean
clean:
	rm -f *.o *~ $(TARGET)
# Command to remove all compiled files to make a clean install