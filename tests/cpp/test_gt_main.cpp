#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>

#include "../../src/PhyloAcc-GT/newick2.h"
#include "../../src/PhyloAcc-common/profile.h"

static void test_load_coal_tree_minimal() {
    std::string tree_path = "tests/data/minimal/tree_coal.tre";
    PhyloTree_theta tree = LoadPhyloTree_theta(tree_path);

    assert(tree.S == 3);
    assert(tree.species_names.size() == 3);
    assert(tree.nodes_names.size() == 5);
    assert(tree.distances.size() == 5);
    assert(std::find(tree.species_names.begin(), tree.species_names.end(), "sp1") != tree.species_names.end());
    assert(std::find(tree.species_names.begin(), tree.species_names.end(), "sp2") != tree.species_names.end());
    assert(std::find(tree.species_names.begin(), tree.species_names.end(), "sp3") != tree.species_names.end());
}

static void test_load_gt_profile_minimal() {
    std::string aln_path = "tests/data/minimal/aln.fa";
    std::string bed_path = "tests/data/minimal/bed.bed";
    PhyloProf prof = LoadPhyloProfiles(aln_path, bed_path);

    assert(prof.S == 3);
    assert(prof.G == 20);
    assert(prof.C == 2);
    assert(prof.element_names.size() == 2);
}

int main() {
    test_load_coal_tree_minimal();
    test_load_gt_profile_minimal();
    std::cout << "GT C++ unit tests passed.\n";
    return 0;
}
