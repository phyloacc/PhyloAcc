#ifndef PHYLOACC_BPP_CONSTRUCTOR_H
#define PHYLOACC_BPP_CONSTRUCTOR_H

#include "profile.h"

#include <armadillo>
#include <set>
#include <string>
#include <vector>

namespace phyloacc
{

struct ElementLayout
{
    std::vector<unsigned int> sizes;
    std::vector<unsigned int> starts;
};

struct SubstitutionEigen
{
    arma::mat eigenvec;
    arma::vec eigenval;
    arma::mat eigeninv;
    arma::vec log_pi;
};

std::vector<int> ParseDelimitedNames(const std::string& names,
                                     const std::vector<std::string>& available_names);

ElementLayout BuildElementLayout(const PhyloProf& profile);

int NumBaseForIndel(double indel);

void ApplyIndelExpansion(double indel, double indel2, arma::mat& submat, arma::vec& pi);

SubstitutionEigen ComputeSubstitutionEigen(const arma::mat& submat, const arma::vec& pi);

void BuildUpperTreeSets(int root,
                        const std::vector<int>& outgroup,
                        const std::vector<int>& conserved_group,
                        int* parent,
                        std::set<int>& upper,
                        std::set<int>& upper_conserve);

std::vector<int> BuildNonOutgroupSubtree(int node_count, const std::set<int>& upper);

}  // namespace phyloacc

#endif
