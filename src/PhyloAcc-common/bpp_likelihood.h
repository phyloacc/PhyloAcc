#ifndef PHYLOACC_BPP_LIKELIHOOD_H
#define PHYLOACC_BPP_LIKELIHOOD_H

#include <armadillo>

#include <cmath>
#include <vector>

namespace phyloacc
{
template <typename LogMultiFn, typename LogExpSumFn>
double ComputeLogLikelihood(std::vector<std::vector<arma::vec>>& lambda,
                            double indel,
                            double indel2,
                            int start_index,
                            int end_index,
                            std::vector<unsigned int>& site_ids,
                            double root_prob,
                            std::vector<int>& subtree,
                            int tip_count,
                            int (*children)[2],
                            double* distances,
                            LogMultiFn log_multi_fn,
                            LogExpSumFn log_exp_sum_fn)
{
    double result = 0;
    arma::mat x(2, 2);

    int root = *subtree.rbegin();
    for (std::vector<int>::iterator it = subtree.begin(); it != subtree.end(); ++it)
    {
        int node = *it;
        if (node < tip_count)
            continue;

        int* node_children = children[node];
        for (int site_index = start_index; site_index < end_index; ++site_index)
            lambda[site_ids[site_index]][node].fill(0);

        for (int child_index = 0; child_index < 2; ++child_index)
        {
            int child = node_children[child_index];
            if (distances[child] > 0)
            {
                double tt = (1 - std::exp(-(indel + indel2) * distances[child])) / (indel + indel2);
                x.at(1, 0) = indel * tt;
                x.at(0, 0) = 1 - x.at(1, 0);
                x.at(0, 1) = indel2 * tt;
                x.at(1, 1) = 1 - x.at(0, 1);
                x = arma::log(x);
            }
            else
            {
                x.fill(-INFINITY);
                x.diag().fill(0);
            }

#pragma omp parallel for schedule(guided)
            for (int site_index = start_index; site_index < end_index; ++site_index)
                lambda[site_ids[site_index]][node] += log_multi_fn(x, lambda[site_ids[site_index]][child]);
        }
    }

    for (int site_index = start_index; site_index < end_index; ++site_index)
    {
        lambda[site_ids[site_index]][root][0] += std::log(1 - root_prob);
        lambda[site_ids[site_index]][root][1] += std::log(root_prob);
        result += log_exp_sum_fn(lambda[site_ids[site_index]][root]);
    }

    return result;
}
}

#endif
