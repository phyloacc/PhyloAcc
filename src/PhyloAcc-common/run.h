#ifndef PHYLOACC_COMMON_RUN_H
#define PHYLOACC_COMMON_RUN_H

#include "config.h"
#include "newick.h"
#include "profile.h"

#include <fstream>
#include <string>
#include <vector>

namespace phyloacc {

enum class ProgramKind {
    ST,
    GT
};

enum class ModelId {
    M0,
    M1,
    M2
};

struct ModelSpec {
    ModelId id;
    int res_z;
    int trace_slot;
    std::string suffix;
};

const ModelSpec& GetModelSpec(ModelId id);
const std::vector<ModelSpec>& AllModelSpecs();

Config LoadConfigForProgram(int argc, char* argv[], ProgramKind kind);
bool ValidateOutputDirectory(const Config& config);
PhyloProf LoadProfile(const Config& config);
PhyloTree LoadSpeciesTree(const Config& config);
void DisplayRunSummary(const PhyloProf& profile, const Config& config, ProgramKind kind);
std::vector<int> ResolveElementIds(const Config& config, int element_count, ProgramKind kind);

struct RunPaths {
    std::string result_folder;
    std::string result_prefix;
    std::string output_prefix;
    std::string output_prefix2;
};

RunPaths MakeRunPaths(const Config& config);

class OutputBundle {
public:
    OutputBundle() = default;
    OutputBundle(const OutputBundle&) = delete;
    OutputBundle& operator=(const OutputBundle&) = delete;
    OutputBundle(OutputBundle&&) = default;
    OutputBundle& operator=(OutputBundle&&) = default;

    std::ofstream hyper;
    std::ofstream likelihood;
    std::ofstream species_names;

    std::ofstream& RatePostZ(ModelId id);
    std::ofstream& Tree(ModelId id);
    void Close();

private:
    std::ofstream rate_m0;
    std::ofstream rate_m1;
    std::ofstream rate_m2;
    std::ofstream tree_m0;
    std::ofstream tree_m1;
    std::ofstream tree_m2;

    friend OutputBundle OpenOutputBundle(ProgramKind kind,
                                         const RunPaths& paths,
                                         const std::vector<std::string>& node_names,
                                         const Config& config);
};

OutputBundle OpenOutputBundle(ProgramKind kind,
                              const RunPaths& paths,
                              const std::vector<std::string>& node_names,
                              const Config& config);

}  // namespace phyloacc

#endif
