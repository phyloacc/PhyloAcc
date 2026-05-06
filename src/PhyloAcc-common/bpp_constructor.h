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

enum class MissingBasePolicy
{
    GapOnly,
    GapNStar
};

struct LeafEncoding
{
    std::vector<std::vector<arma::vec> > lambda;
    std::vector<std::vector<int> > tg;
};

struct ColumnFilterResult
{
    bool filtered;
    int length;
    std::vector<int> removed_sites;
};

struct BppCTraceBuffers
{
    std::vector<double> trace_loglik;
    std::vector<double> trace_full_loglik;
    std::vector<std::vector<int> > trace_z;
    std::vector<double> trace_n_rate;
    std::vector<double> trace_c_rate;
    std::vector<double> trace_l_rate;
    std::vector<double> trace_l2_rate;
    std::vector<double> trace_g_rate;
    std::vector<std::vector<double> > log_emission;
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

int EncodeLeafState(char base, char gapchar, MissingBasePolicy missing_policy);

void EncodeLeafBase(char base,
                    char gapchar,
                    int num_base,
                    MissingBasePolicy missing_policy,
                    arma::vec& lambda,
                    int& tg);

LeafEncoding EncodeLeafAlignment(const std::vector<std::string>& sequences,
                                 int start,
                                 int length,
                                 int species_count,
                                 int node_count,
                                 int num_base,
                                 char gapchar,
                                 MissingBasePolicy missing_policy,
                                 const std::vector<int>& site_order);

std::vector<int> CountMissingBySpecies(const std::vector<std::vector<int> >& tg,
                                       int species_count);

ColumnFilterResult RemoveHighMissingColumns(std::vector<std::vector<arma::vec> >& lambda,
                                            std::vector<std::vector<int> >& tg,
                                            int species_count,
                                            double revgap,
                                            int min_length,
                                            int* simple_block_count = nullptr);

bool IsSimpleOrMissingLeafPattern(std::vector<int> states);

std::vector<bool> BuildMissingNodes(int species_count,
                                    int node_count,
                                    int (*children)[2],
                                    const std::vector<int>& num_missing,
                                    double missing_threshold,
                                    int site_count);

bool ConservedMissingExceeds(const std::vector<bool>& missing,
                             const std::vector<int>& conserved_group,
                             double conserve_prop);

void CollectUpperNodesInSubtree(const std::vector<int>& nodes,
                                const std::set<int>& upper,
                                const std::set<int>& upper_conserve,
                                std::vector<int>& upper_c,
                                std::vector<int>& upper_conserve_c);

BppCTraceBuffers InitializeBppCTraceBuffers(int trace_length,
                                            int node_count,
                                            double initial_l_rate,
                                            double initial_l2_rate,
                                            double initial_g_rate);

}  // namespace phyloacc

#endif
