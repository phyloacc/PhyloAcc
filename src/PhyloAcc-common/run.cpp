#include "run.h"

#include "utils.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

namespace phyloacc {

namespace {

void WriteHyperHeader(std::ofstream& out, const Config& config) {
    out << "iter\tnprior_a\tnprior_b\tcprior_a\tcprior_b\tprior_l_a\tprior_l_b\tprior_g_a\tprior_g_b\n";
    out << 0 << "\t" << config.nprior_a << "\t" << config.nprior_b << "\t"
        << config.cprior_a << "\t" << config.cprior_b << "\t"
        << config.prior_lrate_a << "\t" << config.prior_lrate_b << "\t"
        << config.prior_grate_a << "\t" << config.prior_grate_b << std::endl;
}

void WriteRateHeader(std::ofstream& out, const std::vector<std::string>& node_names) {
    out << "No.\tn_rate\tc_rate\tg_rate\tl_rate\tl2_rate";
    for (std::size_t s = 0; s < node_names.size(); ++s) {
        for (int k = 0; k < 4; ++k) {
            out << "\t" << node_names[s] << "_" << k;
        }
    }
    out << std::endl;
}

std::string ModelPath(const RunPaths& paths, const std::string& stem, ModelId id, const std::string& extension) {
    return paths.output_prefix + stem + GetModelSpec(id).suffix + extension;
}

}  // namespace

const ModelSpec& GetModelSpec(ModelId id) {
    static const ModelSpec m0 = {ModelId::M0, 0, 0, "M0"};
    static const ModelSpec m1 = {ModelId::M1, 2, 1, "M1"};
    static const ModelSpec m2 = {ModelId::M2, 1, 2, "M2"};

    switch (id) {
        case ModelId::M0:
            return m0;
        case ModelId::M1:
            return m1;
        case ModelId::M2:
            return m2;
    }
    return m2;
}

const std::vector<ModelSpec>& AllModelSpecs() {
    static const std::vector<ModelSpec> specs = {
        GetModelSpec(ModelId::M0),
        GetModelSpec(ModelId::M1),
        GetModelSpec(ModelId::M2),
    };
    return specs;
}

Config LoadConfigForProgram(int argc, char* argv[], ProgramKind kind) {
    if (kind == ProgramKind::GT) {
        return LoadConfig(argc, argv, DefaultGTConfig(), true, false);
    }
    return LoadConfig(argc, argv, DefaultSTConfig(), false, true);
}

bool ValidateOutputDirectory(const Config& config) {
    if (!DirectoryExists(config.output_path)) {
        std::cout << "output path doesn't exist or empty!" << std::endl;
        return false;
    }
    return true;
}

PhyloProf LoadProfile(const Config& config) {
    return LoadPhyloProfiles(config.align_path, config.segment_path);
}

PhyloTree LoadSpeciesTree(const Config& config) {
    return LoadPhyloTree(config.phytree_path);
}

void DisplayRunSummary(const PhyloProf& profile, const Config& config, ProgramKind kind) {
    double mean_seg_size = 0;
    for (unsigned int c = 0; c < profile.C; c++) {
        mean_seg_size += (double)(profile.element_pos[c][1] - profile.element_pos[c][0]) / profile.C;
    }
    std::cout << "# total length = " << profile.G << " (" << profile.C << ")"
              << ". # Species = " << profile.S << ". # elements = " << profile.C
              << ". Mean gene set size = " << mean_seg_size << "." << std::endl;
    if (kind == ProgramKind::GT) {
        std::cout << "# Burn-ins = " << config.num_burn * config.num_thin
                  << ". # MCMC Updates = " << config.num_mcmc * config.num_thin
                  << ". # thin = " << config.num_thin << ".  RND SEED = "
                  << config.seed << "." << std::endl;
    } else {
        std::cout << "# Burn-ins = " << config.num_burn
                  << ". # MCMC Updates = " << config.num_mcmc
                  << ". # adaptive frequency = " << config.num_thin
                  << ".  RND SEED = " << config.seed << "." << std::endl;
    }
    std::cout << "# Threads = " << config.num_thread << std::endl << std::endl;
}

std::vector<int> ResolveElementIds(const Config& config, int element_count, ProgramKind kind) {
    std::vector<int> ids;
    if (config.id_path == "") {
        if (config.batch == -1) {
            const int limit = (kind == ProgramKind::GT) ? 500 : element_count;
            for (int c = 0; c < limit; c++) {
                ids.push_back(c);
            }
        } else {
            int temp = std::ceil(element_count / 3);
            for (int c = config.batch * temp; c < (config.batch + 1) * temp; c++) {
                if (c >= element_count) {
                    break;
                }
                ids.push_back(c);
            }
        }
    } else {
        std::ifstream in_params(config.id_path.c_str());
        if (!in_params) {
            std::cerr << "Cannot open the id file: " << config.id_path.c_str() << std::endl;
            std::exit(1);
        }
        std::string line;
        while (std::getline(in_params, line)) {
            std::istringstream line_stream(line);
            std::string tmp;
            line_stream >> tmp;
            tmp = strutils::trim(tmp);
            if (tmp == "") {
                continue;
            }
            ids.push_back(std::atoi(tmp.c_str()));
        }
    }
    return ids;
}

RunPaths MakeRunPaths(const Config& config) {
    RunPaths paths;
    paths.result_folder = config.output_path;
    paths.result_prefix = config.result_prefix;
    paths.output_prefix = config.output_path + "/" + config.result_prefix;
    paths.output_prefix2 = paths.output_prefix;
    return paths;
}

std::ofstream& OutputBundle::RatePostZ(ModelId id) {
    switch (id) {
        case ModelId::M0:
            return rate_m0;
        case ModelId::M1:
            return rate_m1;
        case ModelId::M2:
            return rate_m2;
    }
    return rate_m2;
}

std::ofstream& OutputBundle::Tree(ModelId id) {
    switch (id) {
        case ModelId::M0:
            return tree_m0;
        case ModelId::M1:
            return tree_m1;
        case ModelId::M2:
            return tree_m2;
    }
    return tree_m2;
}

void OutputBundle::Close() {
    rate_m0.close();
    rate_m1.close();
    rate_m2.close();
    tree_m0.close();
    tree_m1.close();
    tree_m2.close();
    hyper.close();
    likelihood.close();
    species_names.close();
}

OutputBundle OpenOutputBundle(ProgramKind kind,
                              const RunPaths& paths,
                              const std::vector<std::string>& node_names,
                              const Config& config) {
    OutputBundle bundle;

    bundle.hyper.open((paths.output_prefix + "_hyper.txt").c_str());
    WriteHyperHeader(bundle.hyper, config);

    if (config.sample_hyper) {
        bundle.likelihood.open((paths.output_prefix + "_elem_lik.txt").c_str());
        bundle.likelihood.precision(8);
        if (kind == ProgramKind::GT) {
            bundle.likelihood << "No.\tID\tloglik_Full\tloglik_Max" << std::endl;
        } else {
            bundle.likelihood << "No.\tID\tloglik_all\tloglik_Max" << std::endl;
        }
    }

    bundle.rate_m0.open(ModelPath(paths, "_rate_postZ_", ModelId::M0, ".txt").c_str());
    bundle.rate_m1.open(ModelPath(paths, "_rate_postZ_", ModelId::M1, ".txt").c_str());
    bundle.rate_m2.open(ModelPath(paths, "_rate_postZ_", ModelId::M2, ".txt").c_str());
    WriteRateHeader(bundle.rate_m0, node_names);
    WriteRateHeader(bundle.rate_m1, node_names);
    WriteRateHeader(bundle.rate_m2, node_names);

    bundle.species_names.open((paths.output_prefix + "_species_names.txt").c_str());
    for (std::size_t s = 0; s < node_names.size(); ++s) {
        bundle.species_names << node_names[s] << std::endl;
    }
    bundle.species_names.close();

    if (kind == ProgramKind::GT) {
        bundle.tree_m0.open(ModelPath(paths, "_tree_", ModelId::M0, ".txt").c_str());
        bundle.tree_m1.open(ModelPath(paths, "_tree_", ModelId::M1, ".txt").c_str());
        bundle.tree_m2.open(ModelPath(paths, "_tree_", ModelId::M2, ".txt").c_str());
        bundle.tree_m0 << "No.\tprop\tgenetree\n";
        bundle.tree_m1 << "No.\tprop\tgenetree\n";
        bundle.tree_m2 << "No.\tprop\tgenetree\n";
    }

    return bundle;
}

}  // namespace phyloacc
