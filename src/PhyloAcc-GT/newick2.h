#ifndef NEWICK2_H_INCLUDED
#define NEWICK2_H_INCLUDED

#include <map>
#include <vector>
#include <string>
#include <armadillo>
//#include "newick.h"

using namespace std;
using namespace arma;

struct PhyloTree_theta
{
    int S;
    vector< string > species_names;     // species names
    vector< string > nodes_names;       // all nodes names  
    vector< double > distances;         // branch distances    
    vector< vector<bool> > dag; 
};


// load the phylogenetic tree
PhyloTree_theta LoadPhyloTree_theta(string params_path);
//void getThetas(PhyloTree tree_sub, PhyloTree tree_coal);
#endif // NEWICK_H_INCLUDED
