#include "bpp_constructor.h"

#include "bpp_tree.h"
#include "utils.h"

#include <algorithm>
#include <cmath>

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

int EncodeLeafState(char base, char gapchar, MissingBasePolicy missing_policy)
{
    switch (base)
    {
        case 'a':
            return 0;
        case 'c':
            return 1;
        case 'g':
            return 2;
        case 't':
            return 3;
        case 'r':
        case 'y':
        case 'k':
        case 'm':
        case 's':
        case 'w':
            return -1;
        default:
            if (base == gapchar || (missing_policy == MissingBasePolicy::GapNStar && (base == 'n' || base == '*')))
            {
                return 4;
            }
            return 5;
    }
}

void EncodeLeafBase(char base,
                    char gapchar,
                    int num_base,
                    MissingBasePolicy missing_policy,
                    arma::vec& lambda,
                    int& tg)
{
    lambda = arma::zeros<arma::vec>(num_base);
    lambda.fill(-INFINITY);
    tg = -1;

    switch (base)
    {
        case 'a':
            lambda[0] = 0;
            tg = 0;
            break;
        case 'c':
            lambda[1] = 0;
            tg = 1;
            break;
        case 'g':
            lambda[2] = 0;
            tg = 2;
            break;
        case 't':
            lambda[3] = 0;
            tg = 3;
            break;
        case 'r':
            lambda[0] = 0;
            lambda[2] = 0;
            break;
        case 'y':
            lambda[1] = 0;
            lambda[3] = 0;
            break;
        case 'k':
            lambda[2] = 0;
            lambda[3] = 0;
            break;
        case 'm':
            lambda[0] = 0;
            lambda[1] = 0;
            break;
        case 's':
            lambda[1] = 0;
            lambda[2] = 0;
            break;
        case 'w':
            lambda[0] = 0;
            lambda[3] = 0;
            break;
        default:
            lambda.fill(0);
            tg = EncodeLeafState(base, gapchar, missing_policy);
            break;
    }
}

LeafEncoding EncodeLeafAlignment(const std::vector<std::string>& sequences,
                                 int start,
                                 int length,
                                 int species_count,
                                 int node_count,
                                 int num_base,
                                 char gapchar,
                                 MissingBasePolicy missing_policy,
                                 const std::vector<int>& site_order)
{
    LeafEncoding encoding;
    encoding.lambda = std::vector<std::vector<arma::vec> >(
        length, std::vector<arma::vec>(node_count, arma::zeros<arma::vec>(num_base)));
    encoding.tg = std::vector<std::vector<int> >(length, std::vector<int>(node_count, -1));

    for (int species = 0; species < species_count; ++species)
    {
        for (int site = 0; site < length; ++site)
        {
            int source_site = site_order.empty() ? site : site_order[site];
            EncodeLeafBase(sequences[species][start + source_site], gapchar, num_base, missing_policy,
                           encoding.lambda[site][species], encoding.tg[site][species]);
        }
    }

    return encoding;
}

std::vector<int> CountMissingBySpecies(const std::vector<std::vector<int> >& tg,
                                       int species_count)
{
    std::vector<int> num_missing(species_count, 0);
    for (int species = 0; species < species_count; ++species)
    {
        for (std::size_t site = 0; site < tg.size(); ++site)
        {
            if (tg[site][species] == 4)
            {
                num_missing[species]++;
            }
        }
    }
    return num_missing;
}

ColumnFilterResult RemoveHighMissingColumns(std::vector<std::vector<arma::vec> >& lambda,
                                            std::vector<std::vector<int> >& tg,
                                            int species_count,
                                            double revgap,
                                            int min_length,
                                            int* simple_block_count)
{
    ColumnFilterResult result;
    result.filtered = false;
    result.length = static_cast<int>(tg.size());

    if (revgap >= 1)
    {
        return result;
    }

    for (int site = 0; site < static_cast<int>(tg.size()); ++site)
    {
        int missing_count = 0;
        for (int species = 0; species < species_count; ++species)
        {
            if (tg[site][species] >= 4)
            {
                missing_count++;
            }
        }
        if (missing_count > species_count * revgap)
        {
            result.removed_sites.push_back(site);
            if (simple_block_count != nullptr && site < *simple_block_count)
            {
                *simple_block_count = *simple_block_count - 1;
            }
        }
    }

    if (static_cast<int>(tg.size() - result.removed_sites.size()) < min_length)
    {
        result.filtered = true;
        return result;
    }

    for (std::vector<int>::reverse_iterator it = result.removed_sites.rbegin();
         it != result.removed_sites.rend(); ++it)
    {
        lambda.erase(lambda.begin() + *it);
        tg.erase(tg.begin() + *it);
    }

    result.length = static_cast<int>(tg.size());
    return result;
}

bool IsSimpleOrMissingLeafPattern(std::vector<int> states)
{
    std::sort(states.begin(), states.end());
    std::vector<int>::iterator unique_end = std::unique(states.begin(), states.end());
    states.resize(std::distance(states.begin(), unique_end));

    if (states.size() == 1)
    {
        return true;
    }

    if (states.size() == 2)
    {
        return states == std::vector<int>{0, 4} || states == std::vector<int>{0, 5} ||
               states == std::vector<int>{1, 4} || states == std::vector<int>{1, 5} ||
               states == std::vector<int>{2, 4} || states == std::vector<int>{2, 5} ||
               states == std::vector<int>{3, 4} || states == std::vector<int>{3, 5} ||
               states == std::vector<int>{4, 5};
    }

    return states == std::vector<int>{0, 4, 5} || states == std::vector<int>{1, 4, 5} ||
           states == std::vector<int>{2, 4, 5} || states == std::vector<int>{3, 4, 5};
}

std::vector<bool> BuildMissingNodes(int species_count,
                                    int node_count,
                                    int (*children)[2],
                                    const std::vector<int>& num_missing,
                                    double missing_threshold,
                                    int site_count)
{
    std::vector<bool> missing(node_count, false);
    for (int node = species_count; node < node_count; ++node)
    {
        int* child_nodes = children[node];
        for (int child_index = 0; child_index < 2; ++child_index)
        {
            int child = child_nodes[child_index];
            if (child < species_count && num_missing[child] > missing_threshold * site_count)
            {
                missing[child] = true;
            }
        }

        if (missing[child_nodes[0]] && missing[child_nodes[1]])
        {
            missing[node] = true;
        }
    }
    return missing;
}

bool ConservedMissingExceeds(const std::vector<bool>& missing,
                             const std::vector<int>& conserved_group,
                             double conserve_prop)
{
    int missing_count = 0;
    for (std::vector<int>::const_iterator it = conserved_group.begin(); it != conserved_group.end(); ++it)
    {
        if (missing[*it])
        {
            missing_count++;
        }
    }
    return missing_count > conserve_prop * conserved_group.size();
}

void CollectUpperNodesInSubtree(const std::vector<int>& nodes,
                                const std::set<int>& upper,
                                const std::set<int>& upper_conserve,
                                std::vector<int>& upper_c,
                                std::vector<int>& upper_conserve_c)
{
    for (std::vector<int>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        if (upper.find(*it) != upper.end())
        {
            upper_c.push_back(*it);
        }

        if (upper_conserve.find(*it) != upper_conserve.end())
        {
            upper_conserve_c.push_back(*it);
        }
    }
}

BppCTraceBuffers InitializeBppCTraceBuffers(int trace_length,
                                            int node_count,
                                            double initial_l_rate,
                                            double initial_l2_rate,
                                            double initial_g_rate)
{
    BppCTraceBuffers traces;
    traces.trace_loglik = std::vector<double>(trace_length, 0);
    traces.trace_full_loglik = std::vector<double>(trace_length, 0);
    traces.trace_z = std::vector<std::vector<int> >(trace_length, std::vector<int>(node_count, 1));
    traces.trace_n_rate = std::vector<double>(trace_length, 0);
    traces.trace_c_rate = std::vector<double>(trace_length, 0);
    traces.trace_l_rate = std::vector<double>(trace_length, initial_l_rate);
    traces.trace_l2_rate = std::vector<double>(trace_length, initial_l2_rate);
    traces.trace_g_rate = std::vector<double>(trace_length, initial_g_rate);
    traces.log_emission = std::vector<std::vector<double> >(node_count, std::vector<double>(3, 0));
    return traces;
}

}  // namespace phyloacc
