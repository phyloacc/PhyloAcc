#include <cstdlib>
#include <iostream>
#include <string>

#include "../../src/PhyloAcc-common/profile.h"

int main() {
    if (std::getenv("EXPECT_FAIL") == nullptr) {
        std::cerr << "EXPECT_FAIL must be set for this helper.\n";
        return 2;
    }
    const char* aln = std::getenv("MALFORMED_ALN");
    const char* bed = std::getenv("MALFORMED_BED");
    if (aln == nullptr || bed == nullptr) {
        std::cerr << "MALFORMED_ALN and MALFORMED_BED must be set.\n";
        return 2;
    }
    LoadPhyloProfiles(std::string(aln), std::string(bed));
    std::cerr << "Expected LoadPhyloProfiles to fail.\n";
    return 3;
}
