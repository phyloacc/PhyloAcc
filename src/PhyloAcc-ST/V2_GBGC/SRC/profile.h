#ifndef PROFILE_H_INCLUDED
#define PROFILE_H_INCLUDED

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>

using namespace std;

struct PhyloProf
{
    unsigned G, S, C; //, P, CG;
    vector< string > species_names;
    vector< string > element_names;
    vector< double* > element_pos;
    vector<string> element_id;
    vector< string> X;
   
};

// load the phylogenetic profile
PhyloProf LoadPhyloProfiles(string profile_path,string segment_path,string segment_ID="");

#endif // PROFILE_H_INCLUDED
