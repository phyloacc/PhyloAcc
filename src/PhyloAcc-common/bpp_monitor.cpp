#include "bpp_monitor.h"

#include <cmath>
#include <gsl/gsl_randist.h>

namespace phyloacc
{
void CopyTraceZ(int m,
                int node_count,
                const std::vector<int>& Z,
                std::vector<std::vector<int> >& trace_Z)
{
    for (int s = 0; s < node_count; ++s)
        trace_Z[m][s] = Z[s];
}

double ComputeBaseFullLogLik(double trace_loglik,
                             double add_loglik,
                             double c_rate,
                             double cprior_a,
                             double cprior_b,
                             int resZ,
                             double n_rate,
                             double nprior_a,
                             double nprior_b)
{
    double full_loglik = trace_loglik + add_loglik;
    full_loglik += log(gsl_ran_gamma_pdf(c_rate, cprior_a, cprior_b));
    if (resZ != 0)
        full_loglik += log(gsl_ran_gamma_pdf(n_rate, nprior_a, nprior_b));
    return full_loglik;
}

bool UpdateMaxState(int m,
                    int burnin,
                    double full_loglik,
                    double& max_loglik,
                    std::vector<int>& max_Z,
                    int& max_m,
                    const std::vector<int>& Z)
{
    if (m < burnin || full_loglik <= max_loglik)
        return false;

    max_loglik = full_loglik;
    max_Z = Z;
    max_m = m;
    return true;
}
}
