#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

#include "../../src/PhyloAcc-common/run.h"

static std::string make_temp_dir() {
    const char* tmpdir = std::getenv("TMPDIR");
    std::string base = (tmpdir == nullptr || std::string(tmpdir).empty()) ? "." : tmpdir;
    std::string templ_string = base + "/phyloacc-run-common-XXXXXX";
    std::vector<char> templ(templ_string.begin(), templ_string.end());
    templ.push_back('\0');
    char* dir = mkdtemp(templ.data());
    assert(dir != nullptr);
    return std::string(dir);
}

static std::string first_line(const std::string& path) {
    std::ifstream in(path.c_str());
    assert(in.good());
    std::string line;
    std::getline(in, line);
    return line;
}

static bool file_exists(const std::string& path) {
    std::ifstream in(path.c_str());
    return in.good();
}

static std::vector<int> ids_from_file(const std::string& path) {
    phyloacc::Config config;
    config.id_path = path;
    return phyloacc::ResolveElementIds(config, 10, phyloacc::ProgramKind::ST);
}

static void test_model_specs() {
    const phyloacc::ModelSpec& m0 = phyloacc::GetModelSpec(phyloacc::ModelId::M0);
    const phyloacc::ModelSpec& m1 = phyloacc::GetModelSpec(phyloacc::ModelId::M1);
    const phyloacc::ModelSpec& m2 = phyloacc::GetModelSpec(phyloacc::ModelId::M2);

    assert(m0.res_z == 0);
    assert(m0.trace_slot == 0);
    assert(m0.suffix == "M0");

    assert(m1.res_z == 2);
    assert(m1.trace_slot == 1);
    assert(m1.suffix == "M1");

    assert(m2.res_z == 1);
    assert(m2.trace_slot == 2);
    assert(m2.suffix == "M2");

    assert(phyloacc::AllModelSpecs().size() == 3);
}

static void test_resolve_element_ids() {
    phyloacc::Config config;

    std::vector<int> st_ids = phyloacc::ResolveElementIds(config, 5, phyloacc::ProgramKind::ST);
    assert((st_ids == std::vector<int>{0, 1, 2, 3, 4}));

    std::vector<int> gt_ids = phyloacc::ResolveElementIds(config, 5, phyloacc::ProgramKind::GT);
    assert(gt_ids.size() == 500);
    assert(gt_ids.front() == 0);
    assert(gt_ids.back() == 499);

    config.batch = 1;
    std::vector<int> batch_ids = phyloacc::ResolveElementIds(config, 9, phyloacc::ProgramKind::ST);
    assert((batch_ids == std::vector<int>{3, 4, 5}));

    const std::string dir = make_temp_dir();
    const std::string id_path = dir + "/ids.txt";
    {
        std::ofstream out(id_path.c_str());
        out << "2\n\n4\n";
    }
    std::vector<int> file_ids = ids_from_file(id_path);
    assert((file_ids == std::vector<int>{2, 4}));
}

static void test_output_headers() {
    const std::vector<std::string> nodes = {"sp1", "sp2"};
    const std::string expected_rate =
        "No.\tn_rate\tc_rate\tg_rate\tl_rate\tl2_rate"
        "\tsp1_0\tsp1_1\tsp1_2\tsp1_3"
        "\tsp2_0\tsp2_1\tsp2_2\tsp2_3";

    phyloacc::Config st_config;
    st_config.output_path = make_temp_dir();
    st_config.result_prefix = "st-run";
    st_config.sample_hyper = true;
    phyloacc::RunPaths st_paths = phyloacc::MakeRunPaths(st_config);
    phyloacc::OutputBundle st_outputs = phyloacc::OpenOutputBundle(
        phyloacc::ProgramKind::ST, st_paths, nodes, st_config);
    st_outputs.Close();

    assert(file_exists(st_paths.output_prefix + "_rate_postZ_M0.txt"));
    assert(file_exists(st_paths.output_prefix + "_rate_postZ_M1.txt"));
    assert(file_exists(st_paths.output_prefix + "_rate_postZ_M2.txt"));
    assert(first_line(st_paths.output_prefix + "_rate_postZ_M0.txt") == expected_rate);
    assert(first_line(st_paths.output_prefix + "_elem_lik.txt") == "No.\tID\tloglik_all\tloglik_Max");
    assert(first_line(st_paths.output_prefix + "_species_names.txt") == "sp1");

    phyloacc::Config gt_config;
    gt_config.output_path = make_temp_dir();
    gt_config.result_prefix = "gt-run";
    gt_config.sample_hyper = true;
    phyloacc::RunPaths gt_paths = phyloacc::MakeRunPaths(gt_config);
    phyloacc::OutputBundle gt_outputs = phyloacc::OpenOutputBundle(
        phyloacc::ProgramKind::GT, gt_paths, nodes, gt_config);
    gt_outputs.Close();

    assert(first_line(gt_paths.output_prefix + "_rate_postZ_M1.txt") == expected_rate);
    assert(first_line(gt_paths.output_prefix + "_tree_M0.txt") == "No.\tprop\tgenetree");
    assert(first_line(gt_paths.output_prefix + "_tree_M1.txt") == "No.\tprop\tgenetree");
    assert(first_line(gt_paths.output_prefix + "_tree_M2.txt") == "No.\tprop\tgenetree");
    assert(first_line(gt_paths.output_prefix + "_elem_lik.txt") == "No.\tID\tloglik_Full\tloglik_Max");
}

int main() {
    test_model_specs();
    test_resolve_element_ids();
    test_output_headers();
    std::cout << "Run common C++ unit tests passed.\n";
    return 0;
}
