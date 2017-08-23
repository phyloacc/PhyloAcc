CXX = g++ -Wall
CXXFLAGS = -g -O2 -std=c++11
LDFLAGS = -lgsl -lm -lgslcblas -larmadillo -fopenmp

GSL_INCLUDE = #/n/sw/fasrcsw/apps/Core/gcc/5.2.0-fasrc02 #/n/sw/fasrcsw/apps/Comp/gcc/4.8.2-fasrc01/gsl/2.1-fasrc01/include #/n/sw/centos6/gsl-1.16_gcc-4.4.7/include
GSL_LIB = #/n/sw/fasrcsw/apps/Comp/gcc/4.8.2-fasrc01/gsl/2.1-fasrc01/bin #/n/sw/centos6/gsl-1.16_gcc-4.4.7/lib

GSL_INCLUDE_MAC = /opt/local/include/
GSL_LIB_MAC = /opt/local/lib/

LSOURCE = main.cpp utils.cpp profile.cpp newick.cpp bpp.cpp bpp_c.cpp bpp_c2.cpp
LHEADER = utils.h newick.h profile.h bpp.hpp bpp_c.hpp

mac: $(LSOURCE) $(HEADER)
	  $(CC) -I$(GSL_INCLUDE_MAC) -L$(GSL_LIB_MAC) $(LSOURCE) -o PhyloAcc $(LDFLAGS)


linux: $(LSOURCE) $(HEADER)
	  $(CXX)   $(CXXFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(LSOURCE) -o PhyloAcc $(LDFLAGS)



clean:
	-rm PhyloAcc
