TARGET=PhyloAcc
# ifeq ($(shell uname),Darwin)
# 	CXX=g++
# else
# 	CXX=g++
# endif
CXX=g++

CFLAGS=-Wall -g -O2 -std=c++11
LDFLAGS=-lgsl -lm -lgslcblas -larmadillo -fopenmp

# FIXME: set variables in shell works. in makefile, cannot recognize the pipline, and only return the first command
# CONDA_ENV_PATH := $(shell conda info -e | grep "\*" | awk -F '[[:space:]]+' '{print $3}')
# CONDA_ENV_INCLUDE := ${CONDA_ENV_PATH}/include/
# CONDA_ENV_LIB := ${CONDA_ENV_PATH}/lib/

SRC_DIR=SRC
SRCS=$(SRC_DIR)/*.cpp
INCLUDES=$(SRC_DIR)/*.h $(SRC_DIR)/*.hpp
PREFIX=/usr/local

$(TARGET): $(SRCS) $(INCLUDES)
	$(CXX) $(CFLAGS) -I${CONDA_ENV_INCLUDE} -L${CONDA_ENV_LIB} $(SRCS) -o $(TARGET) $(LDFLAGS)

# .PHONY: install
# install: $(TARGET)
# 	cp $< $(DESTDIR)$(PREFIX)/bin/$(TARGET)

# .PHONY: uninstall
# uninstall:
# 	rm -f $(DESTDIR)$(PREFIX)/bin/$(TARGET)

.PHONY: clean
clean:
	rm -f *.o *~ $(TARGET)
