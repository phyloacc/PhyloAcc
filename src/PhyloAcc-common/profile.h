#ifndef PROFILE_H
#define PROFILE_H

#include <string>
#include <vector>

using namespace std;

struct PhyloProf
{
    unsigned G, S, C;
    vector<string> species_names;
    vector<string> element_names;
    vector<string> element_tree;
    vector<double*> element_pos;
    vector<string> element_id;
    vector<string> X;
};

PhyloProf LoadPhyloProfiles(string profile_path, string segment_path, string segment_ID = "");

#endif
