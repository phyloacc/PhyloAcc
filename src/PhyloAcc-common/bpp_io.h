#ifndef PHYLOACC_BPP_IO_H
#define PHYLOACC_BPP_IO_H

#include "profile.h"

#include <fstream>
#include <vector>

namespace phyloacc
{
void WriteInitLikelihoods(const PhyloProf& profile,
                          std::ofstream& out_lik,
                          const std::vector<int>& ids,
                          const std::vector<double>& log_liks_sgl,
                          const std::vector<std::vector<double>>& log_liks_Z);

void WriteElemZFiles(const std::string& output_path,
                     int node_count,
                     const std::vector<std::string>& node_names,
                     const std::vector<int>& ids,
                     const std::vector<std::vector<std::vector<int>>>& max_Z,
                     const std::vector<std::vector<std::string>>* genetrees = nullptr);
}

#endif
