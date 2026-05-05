#include "bpp_io.h"

namespace phyloacc
{
void WriteInitLikelihoods(const PhyloProf& profile,
                          std::ofstream& out_lik,
                          const std::vector<int>& ids,
                          const std::vector<double>& log_liks_sgl,
                          const std::vector<std::vector<double>>& log_liks_Z)
{
    for (std::vector<int>::const_iterator it = ids.begin(); it != ids.end(); ++it)
    {
        int cc = *it;
        out_lik << cc << "\t" << profile.element_names[cc] << "\t" << log_liks_sgl[cc] << "\t" << log_liks_Z[1][cc];
        out_lik << std::endl;
    }
}

void WriteElemZFiles(const std::string& output_path,
                     int node_count,
                     const std::vector<std::string>& node_names,
                     const std::vector<int>& ids,
                     const std::vector<std::vector<std::vector<int>>>& max_Z,
                     const std::vector<std::vector<std::string>>* genetrees)
{
    std::ofstream out_z;
    std::string outpath_elem;

    for (int r = 0; r < 3; r++)
    {
        if (r == 2)
            outpath_elem = output_path + "_M" + std::to_string(1) + "_elem_Z.txt";
        else if (r == 1)
            outpath_elem = output_path + "_M" + std::to_string(2) + "_elem_Z.txt";
        else
            outpath_elem = output_path + "_M" + std::to_string(0) + "_elem_Z.txt";

        out_z.open(outpath_elem.c_str());
        out_z << "No.";
        for (int s = 0; s < node_count; s++)
            out_z << "\t" << node_names[s];
        if (genetrees)
            out_z << "\tgenetree";
        out_z << std::endl;

        for (std::vector<int>::const_iterator it = ids.begin(); it != ids.end(); ++it)
        {
            int c = *it;
            out_z << c;
            for (int s = 0; s < node_count; s++)
                out_z << "\t" << max_Z[r][c][s];
            if (genetrees)
                out_z << "\t" << (*genetrees)[r][c];
            out_z << std::endl;
        }
        out_z.close();
    }
}
}
