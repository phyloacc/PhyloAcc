#ifndef PHYLOACC_COMMON_CONFIG_H
#define PHYLOACC_COMMON_CONFIG_H

#include <string>
#include <vector>

namespace phyloacc {

struct Config {
    std::string params_path;
    std::string phytree_path;
    std::string align_path;
    std::string output_path = "";
    std::string output_path2 = "";
    std::string segment_path;
    std::string id_path = "";
    std::string result_prefix = "test";
    std::string tree_coal_unit = "";
    std::string outgroup;
    std::string targetspecies;
    std::string conservegroup;
    std::string deepcoal_species = "";
    double conserve_prop = 0.8;

    int num_thread = 1;
    int num_burn = 200;
    int num_mcmc = 800;
    int num_thin = 500;
    int num_chain = 0;

    double prep_lrate = 0.3;
    double prep_lrate2 = 0.1;
    double prep_grate = 0.5;

    double prior_lrate2_a = 1.0;
    double prior_lrate2_b = 9.0;
    double prior_lrate_a = 1.0;
    double prior_lrate_b = 9.0;
    double prior_grate_a = 3.0;
    double prior_grate_b = 1.0;

    double ratio0 = 0.5;
    double ratio1 = 1.0;
    double missing_thres = 0.8;

    double nprior_a = 10.0;
    double nprior_b = 0.2;
    double cprior_a = 5.0;
    double cprior_b = 0.04;
    int ropt = 1;
    double cub = 1.0;
    double nlb = 0.6;

    int batch = -1;
    int seed = 5;
    int seed2 = 1;
    double indel = 0.0;
    double indel2 = 0.0;
    bool sample_indel = false;
    bool sample_hyper = false;
    char gapchar = '-';
    bool verbose = false;
    double consToMis = 0.01;
    bool prune = false;
    double revgap = 1.0;
    int min_length = 50;

    bool WL = true;
    bool simulate = false;
    bool verboseGT = true;
    int block = 15;
    double br_sample_cutoff = 10.0;
    double theta_cutoff = 1.0;
    std::vector<double> prior_dir_par = std::vector<double>(2, 10.0);
};

Config DefaultSTConfig();
Config DefaultGTConfig();
Config LoadConfig(int argc, char* argv[], Config config, bool is_gt, bool print_reading_message);
bool DirectoryExists(const std::string& path);
bool StringToBool(const std::string& s);

}  // namespace phyloacc

#endif
