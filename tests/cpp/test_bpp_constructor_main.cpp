#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "../../src/PhyloAcc-common/bpp_constructor.h"

static void test_parse_delimited_names() {
    std::vector<std::string> names = {"sp1", "sp2", "anc"};
    std::vector<int> parsed = phyloacc::ParseDelimitedNames("sp2;anc", names);
    assert((parsed == std::vector<int>{1, 2}));
}

static void test_element_layout() {
    double elem0[2] = {10, 15};
    double elem1[2] = {50, 58};

    PhyloProf profile;
    profile.element_names = {"elem0", "elem1"};
    profile.element_pos = {elem0, elem1};

    phyloacc::ElementLayout layout = phyloacc::BuildElementLayout(profile);
    assert((layout.sizes == std::vector<unsigned int>{5, 8}));
    assert((layout.starts == std::vector<unsigned int>{0, 5}));
}

static void test_indel_and_eigen_helpers() {
    assert(phyloacc::NumBaseForIndel(0.0) == 4);
    assert(phyloacc::NumBaseForIndel(0.2) == 5);

    arma::mat submat = {
        {-0.6, 0.2, 0.2, 0.2},
        {0.2, -0.6, 0.2, 0.2},
        {0.2, 0.2, -0.6, 0.2},
        {0.2, 0.2, 0.2, -0.6},
    };
    arma::vec pi(4);
    pi.fill(0.25);

    phyloacc::SubstitutionEigen eigen4 = phyloacc::ComputeSubstitutionEigen(submat, pi);
    assert(eigen4.eigenvec.n_rows == 4);
    assert(eigen4.eigenvec.n_cols == 4);
    assert(eigen4.eigenval.n_elem == 4);
    assert(eigen4.eigeninv.n_rows == 4);
    assert(eigen4.log_pi.n_elem == 4);
    assert(std::fabs(eigen4.log_pi[0] - std::log(0.25)) < 1e-12);

    phyloacc::ApplyIndelExpansion(0.1, 0.05, submat, pi);
    assert(submat.n_rows == 5);
    assert(submat.n_cols == 5);
    assert(pi.n_elem == 5);
    assert(std::fabs(arma::accu(pi) - 1.0) < 1e-10);
    for (std::size_t row = 0; row < submat.n_rows; ++row) {
        assert(std::fabs(arma::accu(submat.row(row))) < 1e-10);
    }
}

static void test_upper_sets_and_subtree() {
    int parent[5] = {3, 3, 4, 4, -1};
    std::set<int> upper;
    std::set<int> upper_conserve;

    phyloacc::BuildUpperTreeSets(4, std::vector<int>{0}, std::vector<int>{1},
                                 parent, upper, upper_conserve);

    assert((upper == std::set<int>{0, 3, 4}));
    assert((upper_conserve == std::set<int>{1, 3, 4}));

    std::vector<int> subtree = phyloacc::BuildNonOutgroupSubtree(5, upper);
    assert((subtree == std::vector<int>{1, 2}));
}

int main() {
    test_parse_delimited_names();
    test_element_layout();
    test_indel_and_eigen_helpers();
    test_upper_sets_and_subtree();
    std::cout << "BPP constructor C++ unit tests passed.\n";
    return 0;
}
