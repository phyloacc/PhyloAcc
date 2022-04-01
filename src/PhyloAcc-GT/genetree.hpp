//  genetree.hpp
//  PhyloAcc_init3-1
//
//  Created by hzr on 2019/5/19.
//  Copyright © 2019 hzr. All rights reserved.
//

#ifndef genetree_h
#define genetree_h

#include "bpp.hpp"

using namespace std;
using namespace arma;

class GTree
{
private:
    int N;
    int S;
    //int seed;
    gsl_rng * RNG;
    
public:
    // gene tree
    int GG;
    int root;
    int    (*children_gene)[2];
    int     *parent_gene;
    bool     *missing_gene;
    double  *heights_gene;
    vector<int> gene_nodes;
    vector<int> var_br_node; 
    vector<int> prob_var_node; //sampling prob of var_br_node. If not uniform, then will be prop to 1/pa(var_br_node)_len; 
    int*  childID_gene; // record the child ID: 0 or 1
    
    vector<vector<int>> temp_coal; // record the coalescents for each species
    
    vector<map<int, int>>  parent_gene2;  // record gene tree plus cross species boundary
    vector<map<int, vector<mat>> > lambda; // N * gene_tree_node * mat
    vector < map<int, vector< int> >> Tg; //current history during each speciation, N*gene_tree_node* GG
    
    
    //GTree(){}
    GTree(int _N, int _S)
    {
        N = _N;
        S = _S;
        
        lambda = vector< map<int, vector<mat>> >(N, map<int, vector<mat>>());
        Tg = vector< map<int, vector < int> >> (N, map<int, vector<int>>());
        
        // gene tree
        children_gene    = new int[N][2];
        parent_gene      = new int[N];
        heights_gene   = new double[N];
        childID_gene = new int[N];
        missing_gene = new bool[N];
        
        parent_gene2 = vector<map<int, int>>(N);
        temp_coal = vector<vector<int>>(N);
    }
    
    GTree(int _N, int _GG, int _S,  gsl_rng* _RNG)
    {
        GG = _GG;
        N = _N;
        S = _S;
        RNG = _RNG;
        
        
        lambda = vector< map<int, vector<mat>> >(N, map<int, vector<mat>>());
        Tg = vector< map<int, vector < int> >> (N, map<int, vector<int>>());
        
        // gene tree
        children_gene    = new int[N][2];
        parent_gene      = new int[N];
        heights_gene   = new double[N];
        childID_gene = new int[N];
        missing_gene = new bool[N];
        
        parent_gene2 = vector<map<int, int>>(N);
        temp_coal = vector<vector<int>>(N);
    }

    ~GTree(){
        delete [] children_gene;
        delete [] parent_gene;
        delete [] heights_gene;
        delete [] childID_gene;
        delete [] missing_gene;
    }
   
        
    void copyto(int len, GTree & gtree);
    void copyfrom(int len, GTree & gtree, int numbase);
    
    //void initTree(string tree_str, vector<bool> & missing, set<int> & upper,BPP & bpp);
    void initTree(string tree_str, BPP & bpp);
    void initTree(vector<bool> & missing, set<int> & upper, BPP & bpp);
    void initTree_Sptop(vector<bool> & missing, set<int> & upper, BPP & bpp);
    void InitTg(int len, BPP & bpp, vector<int> & nodes, int start = 0);
    double priorTree(BPP & bpp);
    void printTree(int s, BPP& bpp, std::stringstream & buffer);
    void printTree2(int s, int & currentS, vector<int> & parent_relabel, vector<double> & branchlen);
    bool Sample_tree(int indicator, BPP & bpp, vector<int> & Z, double n_rate, double c_rate, double consToMis, double nconsToMis, vector<int> lens, vec & loglik,  vec & logp_Z, mat & c_eigenvec, mat & c_eigenval, mat & c_eigeninv, vec& log_pi);

    bool Sample_BranchLen(double delta, int gnode, int indicator, BPP & bpp, vector<int> & Z, double n_rate, double c_rate, double consToMis, double nconsToMis, vector<int> lens, vec & loglik,  vec & logp_Z, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv, double & curGTprior, vec& log_pi);
    void Graft(int branch, int ss, int target, double coal, BPP& bpp, double rate, mat & c_eigenvec, mat & c_eigenval, mat & c_eigeninv);

    int Remove_branch(int branch, int ss, BPP & bpp);
    void getGeneNodes(int s, int & minss);
    double get_logpZ(vector<int> & Z, double consToMis, double nconsToMis);
    double get_logpZ1(BPP & bpp, vector<int> & Z, double consToMis, double nconsToMis);
    
    void Update_Lambda(int start_ss, int branchsib, BPP& bpp, vector<int> & Z, double n_rate, double c_rate, mat & c_eigenvec, mat & c_eigenval, mat & c_eigeninv);
    void Update_Tg(int g, vector<bool> visited, BPP& bpp, bool tosample,vector<int> & nodes, vector<int> & Z, double n_rate, double c_rate, vec log_pi, mat & c_eigenvec, mat & c_eigenval, mat & c_eigeninv, int start =0);
    void Simulate_Tg(int len, BPP& bpp, vector<int> & nodes, vector<int> & Z, double n_rate, double c_rate, vec log_pi, mat & c_eigenvec, mat & c_eigenval, mat & c_eigeninv, int start=0);

    double CompareTg(vector<int> Tg1, vector<int> Tg2, BPP & bpp);

    bool Sample_tree2(int branch, int indicator, BPP & bpp, vector<int> & Z, double n_rate, double c_rate, double consToMis, double nconsToMis, vector<int> lens, vec & loglik,  vec & logp_Z, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv, vec& log_pi);
    //bool Sample_tree2_print(int branch, int indicator, BPP & bpp, vector<int> & Z, double n_rate, double c_rate, double consToMis, double nconsToMis, vector<int> lens, vec & loglik,  vec & logp_Z, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv, vec& log_pi);
    
    void printSptree(BPP& bpp, int l);
    int findSp(int gnode, BPP& bpp);
};

#endif /* genetree_h */
