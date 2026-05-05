#include <cstdlib>
#include <iostream>

#include "../../src/PhyloAcc-GT/newick2.h"

int main() {
    if (std::getenv("EXPECT_FAIL") == nullptr) {
        std::cerr << "EXPECT_FAIL must be set for this helper.\n";
        return 2;
    }
    const char* tree_path = std::getenv("TREE_PATH");
    std::string path = tree_path ? std::string(tree_path) : "tests/data/minimal/does_not_exist.tre";
    LoadPhyloTree_theta(path);
    std::cerr << "Expected LoadPhyloTree_theta to fail.\n";
    return 3;
}
