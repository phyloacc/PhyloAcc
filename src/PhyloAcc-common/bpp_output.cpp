#include "bpp_output.h"

#include <algorithm>
#include <sstream>

namespace phyloacc
{
std::string MakeTracePath(const std::string& output_path2, int resZ, int cc)
{
    return output_path2 + "_mcmc_trace_M" + std::to_string(resZ) + "_" + std::to_string(cc) + ".txt";
}

double MedianInPlace(std::vector<double>& values, std::size_t begin, std::size_t end)
{
    std::sort(values.begin() + begin, values.begin() + end);
    return values[(begin + end) / 2];
}

double MeanRange(const std::vector<double>& values, std::size_t begin, std::size_t end)
{
    double total = 0;
    for (std::size_t i = begin; i < end; ++i)
        total += values[i];
    return total / (end - begin);
}

std::vector<std::vector<int> > CountZStates(int node_count,
                                            std::size_t begin,
                                            std::size_t end,
                                            const std::vector<std::vector<int> >& trace_Z,
                                            const std::vector<bool>& missing,
                                            int missing_count_value)
{
    std::vector<std::vector<int> > countZ(node_count, std::vector<int>(4, 0));
    for (int s = 0; s < node_count; ++s)
    {
        for (std::size_t i = begin; i < end; ++i)
            countZ[s][trace_Z[i][s] + 1]++;

        if (missing[s])
            countZ[s][0] = missing_count_value;
    }
    return countZ;
}

std::string FormatInitSummaryRow(int cc,
                                 double n_rate,
                                 double c_rate,
                                 double g_rate,
                                 double l_rate,
                                 double l2_rate,
                                 const std::vector<std::vector<int> >& countZ,
                                 int denom)
{
    std::ostringstream row;
    row << cc << "\t" << n_rate << "\t" << c_rate << "\t" << g_rate << "\t" << l_rate << "\t" << l2_rate;
    for (std::size_t s = 0; s < countZ.size(); ++s)
    {
        for (int k = 0; k < 4; ++k)
            row << "\t" << static_cast<double>(countZ[s][k]) / denom;
    }
    row << "\n";
    return row.str();
}

void WriteInitSummaryRow(std::ofstream& out_Z,
                         int cc,
                         double n_rate,
                         double c_rate,
                         double g_rate,
                         double l_rate,
                         double l2_rate,
                         const std::vector<std::vector<int> >& countZ,
                         int denom)
{
    out_Z << FormatInitSummaryRow(cc, n_rate, c_rate, g_rate, l_rate, l2_rate, countZ, denom);
}

void WriteNodeHeader(std::ofstream& out,
                     const std::vector<std::string>& node_names)
{
    for (std::size_t s = 0; s < node_names.size(); ++s)
        out << node_names[s] << "\t";
}

void WriteZTraceRow(std::ofstream& out,
                    const std::vector<int>& trace_row)
{
    for (std::size_t s = 0; s < trace_row.size(); ++s)
        out << trace_row[s] << "\t";
}
}
