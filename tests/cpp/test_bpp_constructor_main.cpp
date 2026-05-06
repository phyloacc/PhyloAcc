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

static void test_leaf_base_encoding() {
    arma::vec lambda;
    int tg = -99;

    phyloacc::EncodeLeafBase('a', '-', 5, phyloacc::MissingBasePolicy::GapOnly, lambda, tg);
    assert(tg == 0);
    assert(lambda.n_elem == 5);
    assert(lambda[0] == 0);
    assert(std::isinf(lambda[1]) && lambda[1] < 0);

    phyloacc::EncodeLeafBase('r', '-', 5, phyloacc::MissingBasePolicy::GapOnly, lambda, tg);
    assert(tg == -1);
    assert(lambda[0] == 0);
    assert(lambda[2] == 0);
    assert(std::isinf(lambda[1]) && lambda[1] < 0);

    phyloacc::EncodeLeafBase('-', '-', 5, phyloacc::MissingBasePolicy::GapOnly, lambda, tg);
    assert(tg == 4);
    assert(arma::accu(lambda) == 0);

    phyloacc::EncodeLeafBase('n', '-', 5, phyloacc::MissingBasePolicy::GapOnly, lambda, tg);
    assert(tg == 5);
    assert(arma::accu(lambda) == 0);

    phyloacc::EncodeLeafBase('n', '-', 5, phyloacc::MissingBasePolicy::GapNStar, lambda, tg);
    assert(tg == 4);
    assert(arma::accu(lambda) == 0);

    assert(phyloacc::EncodeLeafState('*', '-', phyloacc::MissingBasePolicy::GapOnly) == 5);
    assert(phyloacc::EncodeLeafState('*', '-', phyloacc::MissingBasePolicy::GapNStar) == 4);
}

static void test_leaf_alignment_encoding_and_missing_counts() {
    std::vector<std::string> sequences = {"acn-", "rg*t"};
    std::vector<int> site_order = {2, 0, 3, 1};

    phyloacc::LeafEncoding encoding = phyloacc::EncodeLeafAlignment(
        sequences, 0, 4, 2, 5, 5, '-', phyloacc::MissingBasePolicy::GapNStar, site_order);

    assert(encoding.lambda.size() == 4);
    assert(encoding.lambda[0].size() == 5);
    assert(encoding.lambda[0][0].n_elem == 5);
    assert(encoding.tg[0][0] == 4);
    assert(encoding.tg[0][1] == 4);
    assert(encoding.tg[1][0] == 0);
    assert(encoding.tg[1][1] == -1);
    assert(encoding.tg[2][0] == 4);
    assert(encoding.tg[2][1] == 3);

    std::vector<int> missing = phyloacc::CountMissingBySpecies(encoding.tg, 2);
    assert((missing == std::vector<int>{2, 1}));
}

static void test_high_missing_column_filter() {
    std::vector<std::string> sequences = {"--ca", "-gca", "t-ca"};
    std::vector<int> identity;
    phyloacc::LeafEncoding encoding = phyloacc::EncodeLeafAlignment(
        sequences, 0, 4, 3, 5, 5, '-', phyloacc::MissingBasePolicy::GapOnly, identity);

    int simple_block_count = 2;
    phyloacc::ColumnFilterResult result = phyloacc::RemoveHighMissingColumns(
        encoding.lambda, encoding.tg, 3, 0.5, 2, &simple_block_count);

    assert(!result.filtered);
    assert(result.length == 2);
    assert((result.removed_sites == std::vector<int>{0, 1}));
    assert(simple_block_count == 1);
    assert(encoding.tg.size() == 2);
    assert(encoding.tg[0][0] == 1);

    phyloacc::LeafEncoding filtered_encoding = phyloacc::EncodeLeafAlignment(
        sequences, 0, 4, 3, 5, 5, '-', phyloacc::MissingBasePolicy::GapOnly, identity);
    phyloacc::ColumnFilterResult filtered = phyloacc::RemoveHighMissingColumns(
        filtered_encoding.lambda, filtered_encoding.tg, 3, 0.5, 3);
    assert(filtered.filtered);
    assert(filtered_encoding.tg.size() == 4);
}

static void test_simple_or_missing_leaf_pattern() {
    assert(phyloacc::IsSimpleOrMissingLeafPattern(std::vector<int>{0, 0, 0}));
    assert(phyloacc::IsSimpleOrMissingLeafPattern(std::vector<int>{0, 4, 5}));
    assert(phyloacc::IsSimpleOrMissingLeafPattern(std::vector<int>{3, 5}));
    assert(!phyloacc::IsSimpleOrMissingLeafPattern(std::vector<int>{0, 1}));
    assert(!phyloacc::IsSimpleOrMissingLeafPattern(std::vector<int>{0, 1, 4}));
}

int main() {
    test_parse_delimited_names();
    test_element_layout();
    test_indel_and_eigen_helpers();
    test_upper_sets_and_subtree();
    test_leaf_base_encoding();
    test_leaf_alignment_encoding_and_missing_counts();
    test_high_missing_column_filter();
    test_simple_or_missing_leaf_pattern();
    std::cout << "BPP constructor C++ unit tests passed.\n";
    return 0;
}
