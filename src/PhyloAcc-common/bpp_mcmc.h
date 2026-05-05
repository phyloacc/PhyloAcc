#ifndef PHYLOACC_BPP_MCMC_H
#define PHYLOACC_BPP_MCMC_H

#include <fstream>
#include <vector>

#include <gsl/gsl_rng.h>

namespace phyloacc
{
void InitializeCommonMCMCStorage(
    int C,
    int N,
    double ratio0,
    double ratio1,
    double ind_lrate,
    double ind_lrate2,
    double ind_grate,
    std::vector<std::vector<std::vector<int>>>& max_Z,
    std::vector<std::vector<std::vector<int>>>& cur_Z,
    std::vector<double>& log_liks_null,
    std::vector<std::vector<double>>& log_liks_Z,
    std::vector<double>& log_liks_sgl,
    std::vector<double>& log_liks_resZ,
    std::vector<double>& log_liks_curZ,
    std::vector<double>& log_liks_propZ,
    std::vector<double>& MH_ratio_gain,
    std::vector<double>& MH_ratio_loss,
    std::vector<double>& cur_crate,
    std::vector<double>& cur_nrate,
    std::vector<double>& cur_lrate,
    std::vector<double>& cur_lrate2,
    std::vector<double>& cur_grate);

void SampleProposal(gsl_rng* rng,
                    int iter,
                    double ind_lrate,
                    double ind_grate,
                    double vlr,
                    double vgr,
                    double nprior_a,
                    double nprior_b,
                    double cprior_a,
                    double cprior_b,
                    double indel,
                    double indel2,
                    double& lrate_prop,
                    double& grate_prop,
                    std::ofstream& output);

void SampleHyperparameters(gsl_rng* rng,
                           int iter,
                           std::vector<int>& ids,
                           std::vector<double>& cur_nrate,
                           std::vector<double>& cur_crate,
                           std::vector<double>& cur_lrate,
                           std::vector<double>& cur_lrate2,
                           std::vector<double>& cur_grate,
                           double& nprior_a,
                           double& nprior_b,
                           double& cprior_a,
                           double& cprior_b,
                           double& prior_l_a,
                           double& prior_l_b,
                           double& prior_l2_a,
                           double& prior_l2_b,
                           double& prior_g_a,
                           double& prior_g_b,
                           std::ofstream& output);
}

#endif
