TARGET=PhyloAcc
ifeq ($(shell uname),Darwin)
	CXX=g++-7
else
	CXX=g++
endif

CFLAGS=-Wall -g -O2 -std=c++11
LDFLAGS=-lgsl -lm -lgslcblas -larmadillo -fopenmp

GSL_HOME=/usr/local/
GSL_INCLUDE=${GSL_HOME}/include/
GSL_LIB=${GSL_HOME}/lib/

SRC_DIR=.
SRCS=$(SRC_DIR)/*.cpp
INCLUDES=$(SRC_DIR)/*.h $(SRC_DIR)/*.hpp
PREFIX=/usr/local

$(TARGET): $(SRCS) $(INCLUDES)
	$(CXX) $(CFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) $(SRCS) -o $(TARGET) $(LDFLAGS)

.PHONY: install
install: $(TARGET)
	cp $< $(DESTDIR)$(PREFIX)/bin/$(TARGET)

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(TARGET)

.PHONY: clean
clean:
	rm -f *.o *~ $(TARGET)
