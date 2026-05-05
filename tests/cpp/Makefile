CXX ?= g++
PREFIX ?= $(CONDA_PREFIX)

CPPFLAGS := -I$(PREFIX)/include
CXXFLAGS := -O2 -g -std=c++14
LDFLAGS := -L$(PREFIX)/lib
LDLIBS := -larmadillo
COMMON_DIR := ../../src/PhyloAcc-common
ST_DIR := ../../src/PhyloAcc-ST
GT_DIR := ../../src/PhyloAcc-GT

ST_BIN := phyloacc_cpp_tests_st
GT_BIN := phyloacc_cpp_tests_gt
FAIL_ST_BIN := phyloacc_cpp_fail_st
FAIL_GT_BIN := phyloacc_cpp_fail_gt
FAIL_ST_PROFILE_BIN := phyloacc_cpp_fail_st_profile
FAIL_GT_PROFILE_BIN := phyloacc_cpp_fail_gt_profile

.PHONY: all clean

all: $(ST_BIN) $(GT_BIN) $(FAIL_ST_BIN) $(FAIL_GT_BIN) $(FAIL_ST_PROFILE_BIN) $(FAIL_GT_PROFILE_BIN)

$(ST_BIN): test_st_main.cpp $(COMMON_DIR)/newick.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(ST_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(GT_BIN): test_gt_main.cpp $(GT_DIR)/newick2.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(GT_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_ST_BIN): fail_st_main.cpp $(COMMON_DIR)/newick.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(ST_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_GT_BIN): fail_gt_main.cpp $(GT_DIR)/newick2.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(GT_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_ST_PROFILE_BIN): fail_st_profile_main.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(ST_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_GT_PROFILE_BIN): fail_gt_profile_main.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(GT_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(ST_BIN) $(GT_BIN) $(FAIL_ST_BIN) $(FAIL_GT_BIN) $(FAIL_ST_PROFILE_BIN) $(FAIL_GT_PROFILE_BIN)
