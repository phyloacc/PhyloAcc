CXX = g++ -Wall
CXXFLAGS = -g -O2 -std=c++11
LDFLAGS = -lgsl -lm -lgslcblas -larmadillo -fopenmp

GSL_INCLUDE = /usr/local/include/
GSL_LIB = /usr/local/lib/

GSL_INCLUDE_MAC = /usr/local/include/
GSL_LIB_MAC = /usr/local/lib/

LSOURCE = main.cpp utils.cpp profile.cpp newick.cpp bpp.cpp bpp_c.cpp bpp_c2.cpp
LHEADER = utils.h newick.h profile.h bpp.hpp bpp_c.hpp

mac: $(LSOURCE) $(HEADER)
	  $(CC) -I$(GSL_INCLUDE_MAC) -L$(GSL_LIB_MAC) $(LSOURCE) -o PhyloAcc $(LDFLAGS)


linux: $(LSOURCE) $(HEADER)
	  $(CXX)   $(CXXFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(LSOURCE) -o PhyloAcc $(LDFLAGS)



clean:
	-rm PhyloAcc
