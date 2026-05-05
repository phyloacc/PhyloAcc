#ifndef PHYLOACC_COMMON_BPP_OUTPUT_H
#define PHYLOACC_COMMON_BPP_OUTPUT_H

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

namespace phyloacc
{
std::string MakeTracePath(const std::string& output_path2, int resZ, int cc);

double MedianInPlace(std::vector<double>& values, std::size_t begin, std::size_t end);
double MeanRange(const std::vector<double>& values, std::size_t begin, std::size_t end);

std::vector<std::vector<int> > CountZStates(int node_count,
                                            std::size_t begin,
                                            std::size_t end,
                                            const std::vector<std::vector<int> >& trace_Z,
                                            const std::vector<bool>& missing,
                                            int missing_count_value);

void WriteInitSummaryRow(std::ofstream& out_Z,
                         int cc,
                         double n_rate,
                         double c_rate,
                         double g_rate,
                         double l_rate,
                         double l2_rate,
                         const std::vector<std::vector<int> >& countZ,
                         int denom);

void WriteNodeHeader(std::ofstream& out,
                     const std::vector<std::string>& node_names);

void WriteZTraceRow(std::ofstream& out,
                    const std::vector<int>& trace_row);
}

#endif
