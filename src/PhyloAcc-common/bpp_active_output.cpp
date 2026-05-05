#include "bpp_active_output.h"

#include "bpp_output.h"

#include <fstream>

namespace phyloacc
{
void WriteSTElemLikelihoodFile(const std::string& output_path,
                               const std::vector<std::string>& element_names,
                               const std::vector<int>& ids,
                               const std::vector<double>& log_liks_null,
                               const std::vector<double>& log_liks_resZ,
                               const std::vector<double>& log_liks_sgl,
                               const std::vector<std::vector<double> >& log_liks_Z)
{
    std::string outpath_elem = output_path + "_elem_lik.txt";
    std::ofstream out_lik(outpath_elem.c_str());
    out_lik.precision(8);
    out_lik << "No.\tID\tloglik_Null\tloglik_Acc\tloglik_Full\tlogBF1\tlogBF2\tloglik_Max_M0\tloglik_Max_M1\tloglik_Max_M2" << std::endl;
    for (std::vector<int>::const_iterator it = ids.begin(); it != ids.end(); ++it)
    {
        const int cc = *it;
        out_lik << cc << "\t" << element_names[cc] << "\t" << log_liks_null[cc] << "\t" << log_liks_resZ[cc] << "\t" << log_liks_sgl[cc] << "\t";
        out_lik << log_liks_resZ[cc] - log_liks_null[cc] << "\t" << log_liks_resZ[cc] - log_liks_sgl[cc];
        out_lik << "\t" << log_liks_Z[0][cc] << "\t" << log_liks_Z[2][cc] << "\t" << log_liks_Z[1][cc] << std::endl;
    }
}

void WriteGTElemLikelihoodFile(const std::string& output_path,
                               const std::vector<std::string>& element_names,
                               const std::vector<int>& ids,
                               const std::vector<std::vector<double> >& log_liks_WL)
{
    std::string outpath_elem = output_path + "_elem_lik.txt";
    std::ofstream out_lik(outpath_elem.c_str());
    out_lik.precision(8);
    out_lik << "No.\tID\tloglik_Null_W\tloglik_Acc_W\tloglik_Full_W\tlogBF1\tlogBF2" << std::endl;
    for (std::vector<int>::const_iterator it = ids.begin(); it != ids.end(); ++it)
    {
        const int cc = *it;
        out_lik << cc << "\t" << element_names[cc] << "\t" << log_liks_WL[0][cc] << "\t" << log_liks_WL[2][cc] << "\t" << log_liks_WL[1][cc] << "\t";
        out_lik << log_liks_WL[2][cc] - log_liks_WL[0][cc] << "\t" << log_liks_WL[2][cc] - log_liks_WL[1][cc] << std::endl;
    }
}

void WriteGTPiModeFiles(const std::string& output_path,
                        const std::vector<int>& ids,
                        const std::vector<std::vector<std::vector<double> > >& cur_pi)
{
    std::ofstream out_pi;
    std::string outpath_elem;
    for (int r = 0; r < 3; r++)
    {
        if (r == 0)
            outpath_elem = output_path + "_M0_Beta_Post_pi_mode.txt";
        else if (r == 2)
            outpath_elem = output_path + "_M1_Beta_Post_pi_mode.txt";
        else
            outpath_elem = output_path + "_M2_Beta_Post_pi_mode.txt";

        out_pi.open(outpath_elem.c_str());
        for (std::vector<int>::const_iterator it = ids.begin(); it != ids.end(); ++it)
        {
            const int c = *it;
            out_pi << c;
            for (int b = 0; b < 4; b++)
                out_pi << "\t" << cur_pi[r][c][b];
            out_pi << std::endl;
        }
        out_pi.close();
    }
}

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
                      const std::vector<std::vector<int> >& trace_Z)
{
    std::string outpath_lik = MakeTracePath(output_path2, resZ, cc);
    std::ofstream out_lik;
    out_lik.precision(8);
    if (iter == 0)
    {
        out_lik.open(outpath_lik.c_str());
        out_lik << "loglik\trate_n\trate_c\tgrate\tlrate\tlrate2\t";
        WriteNodeHeader(out_lik, node_names);
        out_lik << std::endl;
    }
    else
    {
        out_lik.open(outpath_lik.c_str(), std::ios::app);
    }

    for (std::size_t i = 0; i < trace_loglik.size(); i++)
    {
        out_lik << trace_loglik[i] << "\t" << trace_n_rate[i] << "\t" << trace_c_rate[i] << "\t" << trace_g_rate[i] << "\t" << trace_l_rate[i] << "\t" << trace_l2_rate[i] << "\t";
        WriteZTraceRow(out_lik, trace_Z[i]);
        out_lik << std::endl;
    }
}

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
                      const std::vector<std::vector<int> >& trace_Z)
{
    std::string outpath_lik = MakeTracePath(output_path2, resZ, cc);
    std::ofstream out_lik;
    out_lik.precision(8);
    if (iter == 0)
    {
        out_lik.open(outpath_lik.c_str());
        out_lik << "iter\tloglik\tindicator\trate_n\trate_c\tpi_A\tGTtop\t";
        WriteNodeHeader(out_lik, node_names);
        out_lik << "grate\tlrate" << std::endl;
    }
    else
    {
        out_lik.open(outpath_lik.c_str(), std::ios::app);
    }

    for (int i = 0; i < num_mcmc + num_burn; i++)
    {
        out_lik << iter << "\t" << trace_loglik[i] << "\t" << trace_indicator[i] << "\t" << trace_n_rate[i] << "\t" << trace_c_rate[i] << "\t" << trace_pi[i][0] << "\t" << trace_GTtopChg[i] << "\t";
        WriteZTraceRow(out_lik, trace_Z[i]);
        out_lik << trace_g_rate[i] << "\t" << trace_l_rate[i] << "\t" << trace_l2_rate[i] << std::endl;
    }
}
}
