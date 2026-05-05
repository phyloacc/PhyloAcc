#include "config.h"

#include <dirent.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace phyloacc {

namespace {

std::string Trim(const std::string& value) {
    const std::string trim_chars = " \"\t\n";
    const std::size_t first = value.find_first_not_of(trim_chars);
    if (first == std::string::npos) {
        return "";
    }
    const std::size_t last = value.find_last_not_of(trim_chars);
    return value.substr(first, last - first + 1);
}

void ParseCommonParam(const std::string& key, const std::string& value, Config& config, bool& handled) {
    handled = true;

    if (key == "PHYTREE_FILE") {
        config.phytree_path = value;
    } else if (key == "ALIGN_FILE") {
        config.align_path = value;
    } else if (key == "SEG_FILE") {
        config.segment_path = value;
    } else if (key == "ID_FILE") {
        config.id_path = value;
    } else if (key == "BATCH") {
        config.batch = std::stoi(value);
    } else if (key == "RESULT_FOLDER") {
        config.output_path = value;
    } else if (key == "PREFIX") {
        config.result_prefix = value;
    } else if (key == "SEED") {
        config.seed = std::stoi(value);
    } else if (key == "INIT_CONSERVE_RATE") {
        config.ratio0 = std::stod(value);
    } else if (key == "INIT_ACCE_RATE") {
        config.ratio1 = std::stod(value);
    } else if (key == "CONSERVE_PRIOR_A") {
        config.cprior_a = std::stod(value);
    } else if (key == "CONSERVE_PRIOR_B") {
        config.cprior_b = std::stod(value);
    } else if (key == "ACCE_PRIOR_A") {
        config.nprior_a = std::stod(value);
    } else if (key == "ACCE_PRIOR_B") {
        config.nprior_b = std::stod(value);
    } else if (key == "ROPT") {
        config.ropt = std::stoi(value);
    } else if (key == "CUB") {
        config.cub = std::stod(value);
    } else if (key == "NLB") {
        config.nlb = std::stod(value);
    } else if (key == "BURNIN") {
        config.num_burn = std::stoi(value);
    } else if (key == "MCMC") {
        config.num_mcmc = std::stoi(value);
    } else if (key == "INIT_LRATE") {
        config.prep_lrate = std::stod(value);
    } else if (key == "INIT_LRATE2") {
        config.prep_lrate2 = std::stod(value);
    } else if (key == "INIT_GRATE") {
        config.prep_grate = std::stod(value);
    } else if (key == "HYPER_LRATE_A") {
        config.prior_lrate_a = std::stod(value);
    } else if (key == "HYPER_LRATE_B") {
        config.prior_lrate_b = std::stod(value);
    } else if (key == "HYPER_GRATE_A") {
        config.prior_grate_a = std::stod(value);
    } else if (key == "HYPER_GRATE_B") {
        config.prior_grate_b = std::stod(value);
    } else if (key == "HYPER_LRATE2_A") {
        config.prior_lrate2_a = std::stod(value);
    } else if (key == "HYPER_LRATE2_B") {
        config.prior_lrate2_b = std::stod(value);
    } else if (key == "CHAIN") {
        config.num_chain = std::stoi(value);
    } else if (key == "OUTGROUP") {
        config.outgroup = value;
    } else if (key == "TARGETSPECIES") {
        config.targetspecies = value;
    } else if (key == "CONSERVE") {
        config.conservegroup = value;
    } else if (key == "CONSERVE_PROP") {
        config.conserve_prop = std::stod(value);
    } else if (key == "GAP_PROP") {
        config.missing_thres = std::stod(value);
    } else if (key == "CONSTOMIS") {
        config.consToMis = std::stod(value);
    } else if (key == "GAPCHAR") {
        config.gapchar = value[0];
    } else if (key == "PRUNE_TREE") {
        config.prune = StringToBool(value);
    } else if (key == "TRIM_GAP_PERCENT") {
        config.revgap = std::stod(value);
    } else if (key == "MIN_LEN") {
        config.min_length = std::stoi(value);
    } else if (key == "INDEL") {
        config.indel = std::stod(value);
    } else if (key == "INDEL2") {
        config.indel2 = std::stod(value);
    } else if (key == "SAMPLE_INDEL") {
        config.sample_indel = StringToBool(value);
    } else if (key == "SAMPLE_HYPER") {
        config.sample_hyper = StringToBool(value);
    } else if (key == "VERBOSE") {
        config.verbose = StringToBool(value);
    } else if (key == "NUM_THREAD") {
        config.num_thread = std::stoi(value);
    } else {
        handled = false;
    }
}

void ParseSTParam(const std::string& key, const std::string& value, Config& config, bool& handled) {
    ParseCommonParam(key, value, config, handled);
    if (handled) {
        return;
    }

    handled = true;
    if (key == "ADAPT_FREQ") {
        config.num_thin = std::stoi(value);
    } else {
        handled = false;
    }
}

void ParseGTParam(const std::string& key, const std::string& value, Config& config, bool& handled) {
    ParseCommonParam(key, value, config, handled);
    if (handled) {
        return;
    }

    handled = true;
    if (key == "SIMULATE") {
        config.simulate = StringToBool(value);
    } else if (key == "TREE_IN_COALESCENT_UNIT") {
        config.tree_coal_unit = value;
    } else if (key == "SEEDS") {
        config.seed2 = std::stoi(value);
    } else if (key == "THIN") {
        config.num_thin = std::stoi(value);
    } else if (key == "WL") {
        config.WL = StringToBool(value);
    } else if (key == "BLK_WL") {
        config.block = std::stoi(value);
    } else if (key == "BR_SAMPLE_THRESHOLD") {
        config.br_sample_cutoff = std::stod(value);
    } else if (key == "THETA_CUTOFF") {
        config.theta_cutoff = std::stod(value);
    } else if (key == "DEEP_COAL_BRANCH") {
        config.deepcoal_species = value;
    } else if (key == "VERBOSE_GENETREE") {
        config.verboseGT = StringToBool(value);
    } else {
        handled = false;
    }
}

}  // namespace

Config DefaultSTConfig() {
    Config config;
    config.num_thin = 500;
    config.prep_lrate = 0.3;
    config.prep_lrate2 = 0.1;
    config.prep_grate = 0.5;
    config.prior_lrate2_a = 1.0;
    config.prior_lrate2_b = 9.0;
    config.prior_lrate_a = 1.0;
    config.prior_lrate_b = 9.0;
    config.prior_grate_a = 3.0;
    config.prior_grate_b = 1.0;
    config.seed = 5;
    config.consToMis = 0.01;
    config.revgap = 1.0;
    config.prior_dir_par = std::vector<double>(2, 10.0);
    return config;
}

Config DefaultGTConfig() {
    Config config;
    config.num_thin = 1;
    config.prep_lrate = 0.5;
    config.prep_lrate2 = 0.1;
    config.prep_grate = 0.8;
    config.prior_lrate2_a = 1.0;
    config.prior_lrate2_b = 1.0;
    config.prior_lrate_a = 1.0;
    config.prior_lrate_b = 1.0;
    config.prior_grate_a = 1.0;
    config.prior_grate_b = 1.0;
    config.seed = 1;
    config.seed2 = 1;
    config.consToMis = 0.5;
    config.revgap = 0.9;
    config.verboseGT = true;
    config.block = 15;
    config.WL = true;
    config.simulate = false;
    config.br_sample_cutoff = 10.0;
    config.theta_cutoff = 1.0;
    config.prior_dir_par = std::vector<double>(2, 10.0);
    return config;
}

Config LoadConfig(int argc, char* argv[], Config config, bool is_gt, bool print_reading_message) {
    std::cout << "Loading input data and running parameters......" << std::endl;

    if (argc > 1) {
        config.params_path = std::string(argv[1]);
    } else {
        std::cout << "Please specify the path of the parameter file." << std::endl;
        std::exit(1);
    }

    std::cout << "Loading program configurations from " << config.params_path << "......" << std::endl;

    std::ifstream in_params(config.params_path.c_str());
    if (!in_params) {
        std::cerr << "Cannot open the parameters file: " << config.params_path.c_str() << std::endl;
        std::exit(1);
    }

    if (print_reading_message) {
        std::cout << "Reading parameters......" << std::endl;
    }

    std::vector<std::string> param_line;
    std::string word;
    int num_words;
    int line_num = 1;
    std::string line;

    while (std::getline(in_params, line)) {
        param_line.clear();
        std::istringstream line_stream(line);
        while (line_stream >> word) {
            param_line.push_back(word);
        }

        num_words = static_cast<int>(param_line.size());
        if (num_words == 0 || param_line[0][0] == '#') {
            continue;
        }

        if (num_words != 2) {
            std::cerr << std::endl << "Line " << line_num << " in the parameter file is not formatted correctly: " << line << std::endl;
            std::cerr << "Each line should contain a parameter and a value separated by a space and no other whitespace." << std::endl << std::endl;
            std::exit(1);
        }

        bool handled = false;
        if (is_gt) {
            ParseGTParam(param_line[0], param_line[1], config, handled);
        } else {
            ParseSTParam(param_line[0], param_line[1], config, handled);
        }

        if (!handled) {
            std::cout << "Unknown parameter: " << param_line[0] << ", skipping..." << std::endl;
        }

        line_num++;
    }

    if (config.prior_lrate2_a == 0) {
        config.prep_lrate2 = 0;
    }

    config.phytree_path = Trim(config.phytree_path);
    config.align_path = Trim(config.align_path);
    config.output_path = Trim(config.output_path);
    config.segment_path = Trim(config.segment_path);
    config.tree_coal_unit = Trim(config.tree_coal_unit);

    return config;
}

bool DirectoryExists(const std::string& path) {
    if (path.empty()) {
        return false;
    }

    DIR* dir = opendir(path.c_str());
    if (dir == NULL) {
        return false;
    }

    (void)closedir(dir);
    return true;
}

bool StringToBool(const std::string& s) {
    return s == "true" || s == "1";
}

}  // namespace phyloacc
