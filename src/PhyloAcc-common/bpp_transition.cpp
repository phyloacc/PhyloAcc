#include "bpp_transition.h"

#include <cassert>
#include <cmath>

namespace phyloacc
{
void SampleTransitionRates(gsl_rng* rng,
                           const std::vector<int>& nodes,
                           const std::vector<int>& Z,
                           const std::vector<int>& fixZ,
                           const int* parent,
                           int root_sentinel,
                           bool assert_parent_lt_two,
                           double prior_g_a,
                           double prior_g_b,
                           double prior_l_a,
                           double prior_l_b,
                           double prior_l2_a,
                           double prior_l2_b,
                           double& gr,
                           double& lr,
                           double& lr2,
                           std::vector<arma::mat>& log_TM_Int)
{
    arma::mat nZ = arma::zeros(3, 3);
    for (std::vector<int>::const_iterator it = nodes.begin(); it < nodes.end() - 1; ++it)
    {
        int p = parent[*it];
        nZ(Z[*it], Z[p]) += 1;
    }

    gr = gsl_ran_beta(rng, prior_g_a + nZ(1, 0) + nZ(2, 0), prior_g_b + nZ(0, 0));

    for (std::vector<int>::const_iterator it = nodes.begin(); it < nodes.end() - 1; ++it)
    {
        if (fixZ[*it] == 1)
        {
            int p = parent[*it];
            assert(p != root_sentinel);
            if (assert_parent_lt_two)
                assert(Z[p] < 2);
            nZ(Z[*it], Z[p]) -= 1;
        }
    }

    lr = gsl_ran_beta(rng, prior_l_a + nZ(2, 1), prior_l_b + nZ(1, 1));
    if (prior_l2_a == 0)
        lr2 = 0;
    else
        lr2 = gsl_ran_beta(rng, prior_l2_a + nZ(1, 2), prior_l2_b + nZ(2, 2));

    for (std::vector<int>::const_iterator it = nodes.begin(); it < nodes.end(); ++it)
    {
        int s = *it;

        if (fixZ[*it] == 1)
        {
            log_TM_Int[*it](1, 1) = 0;
            log_TM_Int[*it](2, 1) = std::log(0.0);

            log_TM_Int[*it](0, 0) = std::log(1 - gr);
            log_TM_Int[*it](1, 0) = std::log(gr);
            log_TM_Int[*it](2, 0) = std::log(0.0);
        }
        else
        {
            log_TM_Int[s](0, 0) = std::log(1 - gr);
            log_TM_Int[s](1, 0) = std::log(gr);
            log_TM_Int[s](2, 0) = std::log(0.0);

            double y = 1 - lr;
            log_TM_Int[s](1, 1) = std::log(y);
            log_TM_Int[s](2, 1) = std::log(1 - y);

            log_TM_Int[s](1, 2) = std::log(lr2);
            log_TM_Int[s](2, 2) = std::log(1 - lr2);
        }
    }
}
}
