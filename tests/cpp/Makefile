CXX ?= g++
PREFIX ?= $(CONDA_PREFIX)

CPPFLAGS := -I$(PREFIX)/include
CXXFLAGS := -O2 -g -std=c++14
LDFLAGS := -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib
LDLIBS := -larmadillo
COMMON_DIR := ../../src/PhyloAcc-common
ST_DIR := ../../src/PhyloAcc-ST
GT_DIR := ../../src/PhyloAcc-GT

ST_BIN := phyloacc_cpp_tests_st
GT_BIN := phyloacc_cpp_tests_gt
FAIL_ST_BIN := phyloacc_cpp_fail_st
FAIL_GT_BIN := phyloacc_cpp_fail_gt
FAIL_CONFIG_BOOL_BIN := phyloacc_cpp_fail_config_bool
FAIL_ST_PROFILE_BIN := phyloacc_cpp_fail_st_profile
FAIL_GT_PROFILE_BIN := phyloacc_cpp_fail_gt_profile
RUN_COMMON_BIN := phyloacc_cpp_tests_run_common
CONSTRUCTOR_BIN := phyloacc_cpp_tests_bpp_constructor
RNG_BIN := phyloacc_cpp_tests_rng

.PHONY: all clean

all: $(ST_BIN) $(GT_BIN) $(FAIL_ST_BIN) $(FAIL_GT_BIN) $(FAIL_CONFIG_BOOL_BIN) $(FAIL_ST_PROFILE_BIN) $(FAIL_GT_PROFILE_BIN) $(RUN_COMMON_BIN) $(CONSTRUCTOR_BIN) $(RNG_BIN)

$(ST_BIN): test_st_main.cpp $(COMMON_DIR)/newick.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(ST_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(GT_BIN): test_gt_main.cpp $(GT_DIR)/newick2.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(GT_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_ST_BIN): fail_st_main.cpp $(COMMON_DIR)/newick.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(ST_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_GT_BIN): fail_gt_main.cpp $(GT_DIR)/newick2.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(GT_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_CONFIG_BOOL_BIN): fail_config_bool_main.cpp $(COMMON_DIR)/config.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_ST_PROFILE_BIN): fail_st_profile_main.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(ST_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(FAIL_GT_PROFILE_BIN): fail_gt_profile_main.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(GT_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(RUN_COMMON_BIN): test_run_common_main.cpp $(COMMON_DIR)/run.cpp $(COMMON_DIR)/config.cpp $(COMMON_DIR)/newick.cpp $(COMMON_DIR)/profile.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(ST_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(CONSTRUCTOR_BIN): test_bpp_constructor_main.cpp $(COMMON_DIR)/bpp_constructor.cpp $(COMMON_DIR)/bpp_tree.cpp $(COMMON_DIR)/utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) -I$(ST_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(RNG_BIN): test_rng_main.cpp $(COMMON_DIR)/rng.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(COMMON_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(ST_BIN) $(GT_BIN) $(FAIL_ST_BIN) $(FAIL_GT_BIN) $(FAIL_CONFIG_BOOL_BIN) $(FAIL_ST_PROFILE_BIN) $(FAIL_GT_PROFILE_BIN) $(RUN_COMMON_BIN) $(CONSTRUCTOR_BIN) $(RNG_BIN)
