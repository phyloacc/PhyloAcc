#include "bpp_mcmc.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <iostream>

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
    std::vector<double>& cur_grate)
{
    max_Z = std::vector<std::vector<std::vector<int>>>(3, std::vector<std::vector<int>>(C, std::vector<int>(N, 0)));
    cur_Z = std::vector<std::vector<std::vector<int>>>(3, std::vector<std::vector<int>>(C, std::vector<int>(N, 0)));

    log_liks_null = std::vector<double>(C, 0);
    log_liks_Z = std::vector<std::vector<double>>(3, std::vector<double>(C, 0));
    log_liks_sgl = std::vector<double>(C, 0);
    log_liks_resZ = std::vector<double>(C, 0);
    log_liks_curZ = std::vector<double>(C, 0);
    log_liks_propZ = std::vector<double>(C, 0);
    MH_ratio_gain = std::vector<double>(C, 0);
    MH_ratio_loss = std::vector<double>(C, 0);

    cur_crate = std::vector<double>(C, ratio0);
    cur_nrate = std::vector<double>(C, ratio1);
    cur_lrate = std::vector<double>(C, ind_lrate);
    cur_lrate2 = std::vector<double>(C, ind_lrate2);
    cur_grate = std::vector<double>(C, ind_grate);
}
}

void phyloacc::SampleProposal(gsl_rng* rng,
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
                    std::ofstream& output)
{
    lrate_prop = gsl_ran_beta(rng, ind_lrate * vlr, (1 - ind_lrate) * vlr);
    grate_prop = gsl_ran_beta(rng, ind_grate * vgr, (1 - ind_grate - ind_lrate) * vgr);
    grate_prop = grate_prop * (1 - lrate_prop);

    output << iter << "\t" << nprior_a << "\t" << nprior_b << "\t" << cprior_a << "\t" << cprior_b
           << "\t" << indel << "\t" << indel2 << "\t" << ind_grate << "\t" << ind_lrate << std::endl;
}

void phyloacc::SampleHyperparameters(gsl_rng* rng,
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
                           std::ofstream& output)
{
    double p = 1, r = 1;
    double q = 0.1, s = 0.1;

    double vna = 100, vnb = 100, vca = 100, vcb = 100;
    double nprior_a_prop = gsl_ran_gamma(rng, vna, nprior_a / vna);
    double cprior_a_prop = gsl_ran_gamma(rng, vca, cprior_a / vca);

    double nprior_b_prop = gsl_ran_gamma(rng, vnb, nprior_b / vnb);
    double cprior_b_prop = gsl_ran_gamma(rng, vcb, cprior_b / vcb);

    double sum_r = 0;
    double log_prod_r = 0;
    for (std::size_t i = 0; i < ids.size(); i++)
    {
        int c = ids[i];
        sum_r += cur_nrate[c];
        log_prod_r += log(cur_nrate[c]);
    }

    double M_ratio = (nprior_a_prop - 1) * (log(p) + log_prod_r) - (q + sum_r) / nprior_b_prop
                   - (ids.size() + r) * lgamma(nprior_a_prop) - log(nprior_b_prop) * nprior_a_prop * (s + ids.size());
    M_ratio -= (nprior_a - 1) * (log(p) + log_prod_r) - (q + sum_r) / nprior_b
            - (ids.size() + r) * lgamma(nprior_a) - log(nprior_b) * nprior_a * (s + ids.size());

    double H_ratio = log(gsl_ran_gamma_pdf(nprior_a, vna, nprior_a_prop / vna))
                   - log(gsl_ran_gamma_pdf(nprior_a_prop, vna, nprior_a / vna))
                   + log(gsl_ran_gamma_pdf(nprior_b, vnb, nprior_b_prop / vnb))
                   - log(gsl_ran_gamma_pdf(nprior_b_prop, vnb, nprior_b / vnb));

    std::cout << "nrate_MH_ratio: " << M_ratio << ", " << H_ratio << ", " << nprior_a << ", " << nprior_a_prop
              << ", " << nprior_b << ", " << nprior_b_prop << std::endl;

    if (log(gsl_rng_uniform(rng)) < M_ratio + H_ratio)
    {
        nprior_a = nprior_a_prop;
        nprior_b = nprior_b_prop;
    }

    sum_r = 0;
    log_prod_r = 0;
    for (std::size_t i = 0; i < ids.size(); i++)
    {
        int c = ids[i];
        sum_r += cur_crate[c];
        log_prod_r += log(cur_crate[c]);
    }

    M_ratio = (cprior_a_prop - 1) * (log(p) + log_prod_r) - (q + sum_r) / cprior_b_prop
            - (ids.size() + r) * lgamma(cprior_a_prop) - log(cprior_b_prop) * cprior_a_prop * (s + ids.size())
            - ((cprior_a - 1) * (log(p) + log_prod_r) - (q + sum_r) / cprior_b
            - (ids.size() + r) * lgamma(cprior_a) - log(cprior_b) * cprior_a * (s + ids.size()));

    H_ratio = log(gsl_ran_gamma_pdf(cprior_a, vca, cprior_a_prop / vca))
            - log(gsl_ran_gamma_pdf(cprior_a_prop, vca, cprior_a / vca))
            + log(gsl_ran_gamma_pdf(cprior_b, vcb, cprior_b_prop / vcb))
            - log(gsl_ran_gamma_pdf(cprior_b_prop, vcb, cprior_b / vcb));

    std::cout << "crate_MH_ratio: " << M_ratio + H_ratio << ", " << cprior_a << ", " << cprior_a_prop
              << ", " << cprior_b << ", " << cprior_b_prop << std::endl;

    if (log(gsl_rng_uniform(rng)) < M_ratio + H_ratio)
    {
        cprior_a = cprior_a_prop;
        cprior_b = cprior_b_prop;
    }

    double u = gsl_rng_uniform(rng) * (1.3 - 0.7) + 0.7;
    double prior_l_a_prop = prior_l_a * u;

    u = gsl_rng_uniform(rng) * (1.3 - 0.7) + 0.7;
    double prior_l_b_prop = prior_l_b * u;
    double log_p = 0, log_pc = 0;
    for (std::size_t i = 0; i < ids.size(); i++)
    {
        int c = ids[i];
        log_p += log(cur_lrate[c]);
        log_pc += log(1 - cur_lrate[c]);
    }

    M_ratio = (prior_l_a_prop - prior_l_a) * (log_p - 1) + (prior_l_b_prop - prior_l_b) * (log_pc - 1);
    M_ratio += ids.size() * (gsl_sf_lnbeta(prior_l_a, prior_l_b) - gsl_sf_lnbeta(prior_l_a_prop, prior_l_b_prop));

    H_ratio = log(prior_l_a) - log(prior_l_a_prop) + log(prior_l_b) - log(prior_l_b_prop);

    std::cout << "lrate_MH_ratio: " << M_ratio + H_ratio << ", " << prior_l_a << ", " << prior_l_a_prop
              << ", " << prior_l_b << ", " << prior_l_b_prop << std::endl;

    if (log(gsl_rng_uniform(rng)) < M_ratio + H_ratio)
    {
        prior_l_a = prior_l_a_prop;
        prior_l_b = prior_l_b_prop;
    }

    u = gsl_rng_uniform(rng) * (1.3 - 0.7) + 0.7;
    double prior_l2_a_prop = prior_l2_a * u;

    u = gsl_rng_uniform(rng) * (1.3 - 0.7) + 0.7;
    double prior_l2_b_prop = prior_l2_b * u;
    log_p = 0;
    log_pc = 0;
    for (std::size_t i = 0; i < ids.size(); i++)
    {
        int c = ids[i];
        log_p += log(cur_lrate2[c]);
        log_pc += log(1 - cur_lrate2[c]);
    }

    M_ratio = (prior_l2_a_prop - prior_l2_a) * (log_p - 1) + (prior_l2_b_prop - prior_l2_b) * (log_pc - 1);
    M_ratio += ids.size() * (gsl_sf_lnbeta(prior_l2_a, prior_l2_b) - gsl_sf_lnbeta(prior_l2_a_prop, prior_l2_b_prop));

    H_ratio = log(prior_l2_a) - log(prior_l2_a_prop) + log(prior_l2_b) - log(prior_l2_b_prop);

    std::cout << "lrate2_MH_ratio: " << M_ratio + H_ratio << ", " << prior_l2_a << ", " << prior_l2_a_prop
              << ", " << prior_l2_b << ", " << prior_l2_b_prop << std::endl;

    if (log(gsl_rng_uniform(rng)) < M_ratio + H_ratio)
    {
        prior_l2_a = prior_l2_a_prop;
        prior_l2_b = prior_l2_b_prop;
    }

    u = gsl_rng_uniform(rng) * (1.3 - 0.7) + 0.7;
    double prior_g_a_prop = prior_g_a * u;

    u = gsl_rng_uniform(rng) * (1.3 - 0.7) + 0.7;
    double prior_g_b_prop = prior_g_b * u;
    log_p = 0;
    log_pc = 0;
    for (std::size_t i = 0; i < ids.size(); i++)
    {
        int c = ids[i];
        log_p += log(cur_grate[c]);
        log_pc += log(1 - cur_grate[c]);
    }

    M_ratio = (prior_g_a_prop - prior_g_a) * (log_p - 1) + (prior_g_b_prop - prior_g_b) * (log_pc - 1);
    M_ratio += ids.size() * (gsl_sf_lnbeta(prior_g_a, prior_g_b) - gsl_sf_lnbeta(prior_g_a_prop, prior_g_b_prop));

    H_ratio = log(prior_g_a) - log(prior_g_a_prop) + log(prior_g_b) - log(prior_g_b_prop);

    std::cout << "grate_MH_ratio: " << M_ratio + H_ratio << ", " << prior_g_a << ", " << prior_g_a_prop
              << ", " << prior_g_b << ", " << prior_g_b_prop << std::endl;

    if (log(gsl_rng_uniform(rng)) < M_ratio + H_ratio)
    {
        prior_g_a = prior_g_a_prop;
        prior_g_b = prior_g_b_prop;
    }

    output << iter << "\t" << nprior_a << "\t" << nprior_b << "\t" << cprior_a << "\t" << cprior_b
           << "\t" << prior_l_a << "\t" << prior_l_b << "\t" << prior_g_a << "\t" << prior_g_b
           << "\t" << prior_l2_a << "\t" << prior_l2_b << std::endl;
}
