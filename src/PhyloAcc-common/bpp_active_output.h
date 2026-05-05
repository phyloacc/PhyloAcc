#ifndef PHYLOACC_COMMON_BPP_ACTIVE_OUTPUT_H
#define PHYLOACC_COMMON_BPP_ACTIVE_OUTPUT_H

#include <fstream>
#include <string>
#include <vector>

namespace phyloacc
{
void WriteSTElemLikelihoodFile(const std::string& output_path,
                               const std::vector<std::string>& element_names,
                               const std::vector<int>& ids,
                               const std::vector<double>& log_liks_null,
                               const std::vector<double>& log_liks_resZ,
                               const std::vector<double>& log_liks_sgl,
                               const std::vector<std::vector<double> >& log_liks_Z);

void WriteGTElemLikelihoodFile(const std::string& output_path,
                               const std::vector<std::string>& element_names,
                               const std::vector<int>& ids,
                               const std::vector<std::vector<double> >& log_liks_WL);

void WriteGTPiModeFiles(const std::string& output_path,
                        const std::vector<int>& ids,
                        const std::vector<std::vector<std::vector<double> > >& cur_pi);

void WriteSTTraceFile(int iter,
                      const std::string& output_path2,
                      int resZ,
                      int cc,
                      const std::vector<std::string>& node_names,
                      const std::vector<double>& trace_loglik,
                      const std::vector<double>& trace_n_rate,
                      const std::vector<double>& trace_c_rate,
                      const std::vector<double>& trace_g_rate,
                      const std::vector<double>& trace_l_rate,
                      const std::vector<double>& trace_l2_rate,
                      const std::vector<std::vector<int> >& trace_Z);

void WriteGTTraceFile(int iter,
                      const std::string& output_path2,
                      int resZ,
                      int cc,
                      int num_burn,
                      int num_mcmc,
                      const std::vector<std::string>& node_names,
                      const std::vector<double>& trace_loglik,
                      const std::vector<double>& trace_indicator,
                      const std::vector<double>& trace_n_rate,
                      const std::vector<double>& trace_c_rate,
                      const std::vector<std::vector<double> >& trace_pi,
                      const std::vector<int>& trace_GTtopChg,
                      const std::vector<double>& trace_g_rate,
                      const std::vector<double>& trace_l_rate,
                      const std::vector<double>& trace_l2_rate,
                      const std::vector<std::vector<int> >& trace_Z);
}

#endif
