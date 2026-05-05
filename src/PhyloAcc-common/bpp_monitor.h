#ifndef PHYLOACC_COMMON_BPP_MONITOR_H
#define PHYLOACC_COMMON_BPP_MONITOR_H

#include <vector>

namespace phyloacc
{
void CopyTraceZ(int m,
                int node_count,
                const std::vector<int>& Z,
                std::vector<std::vector<int> >& trace_Z);

double ComputeBaseFullLogLik(double trace_loglik,
                             double add_loglik,
                             double c_rate,
                             double cprior_a,
                             double cprior_b,
                             int resZ,
                             double n_rate,
                             double nprior_a,
                             double nprior_b);

bool UpdateMaxState(int m,
                    int burnin,
                    double full_loglik,
                    double& max_loglik,
                    std::vector<int>& max_Z,
                    int& max_m,
                    const std::vector<int>& Z);
}

#endif
