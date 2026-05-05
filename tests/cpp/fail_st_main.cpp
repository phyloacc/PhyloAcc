#include <cstdlib>
#include <iostream>
#include <string>

#include "../../src/PhyloAcc-common/newick.h"

int main() {
    if (std::getenv("EXPECT_FAIL") == nullptr) {
        std::cerr << "EXPECT_FAIL must be set for this helper.\n";
        return 2;
    }
    const char* tree_path = std::getenv("TREE_PATH");
    std::string path = tree_path ? std::string(tree_path) : "tests/data/minimal/does_not_exist.mod";
    LoadPhyloTree(path);
    std::cerr << "Expected LoadPhyloTree to fail.\n";
    return 3;
}
