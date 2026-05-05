#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include "../../src/PhyloAcc-ST/newick.h"
#include "../../src/PhyloAcc-ST/profile.h"

static void test_load_tree_minimal() {
    std::string mod_path = "tests/data/minimal/model.mod";
    PhyloTree tree = LoadPhyloTree(mod_path);

    assert(tree.S == 3);
    assert(tree.species_names.size() == 3);
    assert(tree.nodes_names.size() == 5);
    assert(tree.distances.size() == 5);
    assert(tree.subs_rate.n_rows == 4);
    assert(tree.subs_rate.n_cols == 4);

    // Ensure rate matrix has non-zero diagonal magnitudes.
    for (int i = 0; i < 4; i++) {
        assert(std::fabs(tree.subs_rate(i, i)) > 0.0);
    }
}

static void test_load_profile_minimal() {
    std::string aln_path = "tests/data/minimal/aln.fa";
    std::string bed_path = "tests/data/minimal/bed.bed";
    PhyloProf prof = LoadPhyloProfiles(aln_path, bed_path);

    assert(prof.S == 3);
    assert(prof.G == 20);
    assert(prof.C == 2);
    assert(prof.element_pos.size() == 2);
    assert(static_cast<int>(prof.element_pos[0][0]) == 0);
    assert(static_cast<int>(prof.element_pos[0][1]) == 10);
    assert(static_cast<int>(prof.element_pos[1][0]) == 10);
    assert(static_cast<int>(prof.element_pos[1][1]) == 20);
}

int main() {
    test_load_tree_minimal();
    test_load_profile_minimal();
    std::cout << "C++ unit tests passed.\n";
    return 0;
}
