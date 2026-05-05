#include "bpp_data.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace phyloacc
{
void MatchProfileToTree(PhyloProf& profile, PhyloTree& tree)
{
    bool success_match = true;
    int species_count = tree.S;
    int profile_species_count = profile.S;
    vector<int> reorder(species_count);

    for (int tree_index = 0; tree_index < species_count; tree_index++)
    {
        bool has_same_species = false;
        string tree_name = tree.species_names[tree_index];
        for (int profile_index = 0; profile_index < profile_species_count; profile_index++)
        {
            string profile_name = profile.species_names[profile_index];
            if (tree_name == profile_name)
            {
                has_same_species = true;
                reorder[tree_index] = profile_index;
                break;
            }
        }
        if (!has_same_species)
        {
            cout << "No matrix species " << profile.species_names[tree_index] << " found in tree." << endl;
            success_match = false;
            break;
        }
    }

    if (!success_match)
    {
        cout << endl << "The species in phylogenetic profile and tree cannot be matched literally:" << endl;
        cout << "The program will use the default mapping in data:" << endl;
        for (int species_index = 0; species_index < species_count; species_index++)
            cout << "(" << profile.species_names[species_index] << "\t=  " << tree.species_names[species_index] << ")" << endl;
        cout << endl;
        return;
    }

    cout << "The species in profile and tree match perfectly. Reorder the species in profile matrix by the tree." << endl << endl;
    vector<string> old_X = profile.X;
    for (int species_index = 0; species_index < species_count; species_index++)
    {
        int reorder_species_index = reorder[species_index];
        profile.X[species_index] = old_X[reorder_species_index];
    }
    profile.species_names = tree.species_names;
}
}
