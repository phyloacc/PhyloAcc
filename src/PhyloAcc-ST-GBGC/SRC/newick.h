#ifndef NEWICK_H_INCLUDED
#define NEWICK_H_INCLUDED

#include <vector>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

struct PhyloTree
{
    int S;                              // num of living species
    vector< string > species_names;     // species names
    vector< string > nodes_names;       // all nodes names
    vector< vector<bool> > dag;         // connection matrix
    vector< double > distances;         // branch distances
    // base compostion and substritution rate
    vec pi;
    mat subs_rate;
};

// load the phylogenetic tree
PhyloTree LoadPhyloTree(string params_path);

#endif // NEWICK_H_INCLUDED
