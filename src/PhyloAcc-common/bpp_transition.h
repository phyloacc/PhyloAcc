#ifndef PHYLOACC_COMMON_BPP_TRANSITION_H
#define PHYLOACC_COMMON_BPP_TRANSITION_H

#include <armadillo>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <vector>

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
                           std::vector<arma::mat>& log_TM_Int);
}

#endif
