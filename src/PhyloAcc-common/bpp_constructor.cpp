#include "bpp_constructor.h"

#include "bpp_tree.h"
#include "utils.h"

#include <algorithm>

namespace phyloacc
{

std::vector<int> ParseDelimitedNames(const std::string& names,
                                     const std::vector<std::string>& available_names)
{
    std::vector<int> ids;
    std::vector<std::string> items = strutils::split(strutils::trim(names), ';');
    for (std::vector<std::string>::const_iterator it = items.begin(); it != items.end(); ++it)
    {
        std::ptrdiff_t pos = std::find(available_names.begin(), available_names.end(), *it) - available_names.begin();
        ids.push_back(static_cast<int>(pos));
    }
    return ids;
}

ElementLayout BuildElementLayout(const PhyloProf& profile)
{
    ElementLayout layout;
    layout.sizes = std::vector<unsigned int>(profile.element_names.size());
    layout.starts = std::vector<unsigned int>(profile.element_names.size());

    int start = 0;
    for (std::size_t c = 0; c < profile.element_names.size(); ++c)
    {
        int end = start + profile.element_pos[c][1] - profile.element_pos[c][0];
        int len = end - start;
        layout.sizes[c] = static_cast<unsigned int>(len);
        layout.starts[c] = static_cast<unsigned int>(start);
        start = end;
    }
    return layout;
}

int NumBaseForIndel(double indel)
{
    if (indel < 1e-10)
    {
        return 4;
    }
    return 5;
}

void ApplyIndelExpansion(double indel, double indel2, arma::mat& submat, arma::vec& pi)
{
    submat *= (1 - indel);
    arma::mat insertion = arma::ones<arma::mat>(4, 1) * indel;
    arma::mat deletion = arma::ones<arma::mat>(1, 5) * indel2;

    submat.insert_cols(4, insertion);
    submat.insert_rows(4, deletion);

    arma::colvec row_sums = arma::sum(submat, 1);
    submat.diag() -= row_sums;
    arma::mat stationary = arma::null(submat.t());
    pi = stationary / arma::accu(stationary);
}

SubstitutionEigen ComputeSubstitutionEigen(const arma::mat& submat, const arma::vec& pi)
{
    arma::cx_mat bvec;
    arma::cx_vec aval;
    arma::eig_gen(aval, bvec, submat);

    SubstitutionEigen eigen;
    eigen.eigenval = arma::conv_to<arma::vec>::from(aval);
    eigen.eigenvec = arma::conv_to<arma::mat>::from(bvec).t();
    eigen.eigeninv = arma::inv(eigen.eigenvec);
    eigen.log_pi = arma::log(pi);
    return eigen;
}

void BuildUpperTreeSets(int root,
                        const std::vector<int>& outgroup,
                        const std::vector<int>& conserved_group,
                        int* parent,
                        std::set<int>& upper,
                        std::set<int>& upper_conserve)
{
    CollectUpperTreeNodes(root, outgroup, parent, upper);
    CollectUpperTreeNodes(root, conserved_group, parent, upper_conserve);
}

std::vector<int> BuildNonOutgroupSubtree(int node_count, const std::set<int>& upper)
{
    std::vector<int> subtree;
    for (int node = 0; node < node_count - 1; ++node)
    {
        if (upper.find(node) == upper.end())
        {
            subtree.push_back(node);
        }
    }
    return subtree;
}

}  // namespace phyloacc
