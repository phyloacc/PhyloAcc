TARGET=PhyloAcc-GT_piQ
CFLAGS=-Wall -g -O2 -std=c++11
LDFLAGS=-lgsl -lm -lgslcblas -larmadillo -fopenmp

GSL_INCLUDE=${GSL_HOME}/include
GSL_LIB=${GSL_HOME}/lib64

SRC_DIR=.
SRCS=$(SRC_DIR)/*.cpp
INCLUDES=$(SRC_DIR)/*.h $(SRC_DIR)/*.hpp
PREFIX=/usr/local
ARMA900_INCLUDE=/n/home05/hyan/armadillo-9.900.1/include
ARMA900_LIB=/n/home05/hyan/armadillo-9.900.1

$(TARGET): $(SRCS) $(INCLUDES)
	$(CXX) $(CFLAGS) -I$(GSL_INCLUDE) -I$(ARMA900_INCLUDE) -L$(GSL_LIB) -L$(ARMA900_LIB) $(SRCS) -o $(TARGET) $(LDFLAGS)

.PHONY: install
install: $(TARGET)
	cp $< $(DESTDIR)$(PREFIX)/bin/$(TARGET)

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(TARGET)

.PHONY: clean
clean:
	rm -f *.o *~ $(TARGET)
