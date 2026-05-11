#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

#include "../../src/PhyloAcc-common/config.h"

static std::string make_temp_dir() {
    const char* tmpdir = std::getenv("TMPDIR");
    std::string base = (tmpdir == nullptr || std::string(tmpdir).empty()) ? "." : tmpdir;
    std::string templ_string = base + "/phyloacc-config-bool-XXXXXX";
    std::vector<char> templ(templ_string.begin(), templ_string.end());
    templ.push_back('\0');
    char* dir = mkdtemp(templ.data());
    if (dir == nullptr) {
        std::cerr << "Could not create temporary directory.\n";
        std::exit(2);
    }
    return std::string(dir);
}

int main() {
    if (std::getenv("EXPECT_FAIL") == nullptr) {
        std::cerr << "EXPECT_FAIL must be set for this helper.\n";
        return 2;
    }

    const std::string dir = make_temp_dir();
    const std::string cfg_path = dir + "/bad-bool.cfg";
    {
        std::ofstream out(cfg_path.c_str());
        out << "WL maybe\n";
    }

    std::vector<char> arg0;
    arg0.push_back('t');
    arg0.push_back('\0');
    std::vector<char> arg1(cfg_path.begin(), cfg_path.end());
    arg1.push_back('\0');
    char* argv[] = {arg0.data(), arg1.data()};
    phyloacc::LoadConfig(2, argv, phyloacc::DefaultGTConfig(), true, false);

    std::cerr << "Expected LoadConfig to fail on invalid boolean value.\n";
    return 3;
}
