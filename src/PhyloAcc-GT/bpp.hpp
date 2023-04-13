//  bpp.hpp
//  PhyloAcc
//
//  Created by hzr on 3/8/16.
//  Copyright © 2016 hzr. All rights reserved.
//

#ifndef bpp_hpp
#define bpp_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
#include <armadillo>
#include <queue>

#include <cmath>
#include <sys/time.h>

#include <map>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>

#include "newick.h"
#include "profile.h"
#include "utils.h"

using namespace std;
using namespace arma;


const double LOG_ZERO = -INFINITY;

class BPP
{
    friend class BPP_C;

private:

    int G;  //total base pairs
    int S;
    vector <unsigned int> element_size;
    vector <unsigned int> element_start;
    vector<string> element_tree;

    int root_ingrp;

    double ratio0;
    double ratio1;
    double cub;
    double nlb;
    int ropt;

    // names of species
    vector<int> target_species;
    vector<int> conservedgroup;
    vector<int> dp_species;
    bool cons_sample;
    //int refspecies;

    double conserve_prop;

    vector<int> outgroup;
    set<int> upper;
    set<int> upper_conserve;

    vector<int> subtree;
    vector<int> tosample;

    // MCMC running parameters
    int num_burn;   // num of burn-in updates
    int num_mcmc;   // num of MCMC updates
    int num_thin;   // num of updates between two samples


    vector < vector< string > > X;  //element * current profile

    vec pi; //used vec since newick uses vec
    vec log_pi;
    vector<double> prior_dir_param;
    vec inst_rate; //instantaneous rate of substitution

    mat eigenvec; vec eigenval;
    mat eigeninv;
    mat submat;
    double indel;
    double indel2;

    mat eigenvecprop; vec eigenvalprop;
    mat eigeninvprop; vec piprop;

    int moveroot = -10; // move root to one of its children

    vector<vector < vector <int > >> Max_Z;
    vector<vector < vector <int > >> cur_Z;
    
    vector<vector<string>> genetrees;

    double ind_lrate, ind_lrate2;  //for Z
    double ind_grate;
    double vlr = 100; // proposal variance of loss rate
    double vgr = 100; // proposal variance of gain rate
    //double prior_glr[3];
    double prior_l_a, prior_l_b;
    double prior_l2_a, prior_l2_b;
    double prior_g_a, prior_g_b;

    vector<mat >  TM_Int;
    vector<mat> log_TM_Int;

    vector<mat > log_cache_TM_null;

    vector< double >  log_liks_sgl;
    vector< double >  log_liks_sgl_L;
   
    vector< double > log_liks_curZ; //Store current loglik
    vector< double > log_liks_propZ; //Store loglik with prop indel parameter

    vector< double > MH_ratio_gain; //
    vector< double > MH_ratio_loss; //

    vector< double > log_liks_null;
    vector< double > log_liks_resZ;
    vector< double >  log_liks_null_L;
    vector< double >  log_liks_resZ_L;

    vector< double > cur_nrate;
    vector< double > cur_crate;

    vector< double > cur_lrate;
    vector< double > cur_grate;
    vector< double > cur_lrate2;

    //Han*: save running max pi under MCMC; also acgt counts
    vector<vector<vector<double>>> cur_pi; // resZ=3 * C * pi

    vector < vector <mat>> log_cache_TM_neut_ensemble;
    vector < vector <mat>> log_cache_TM_cons_ensemble;
    vector<double> neut_prior_density;
    vector<double> cons_prior_density;

    //hyperparameter of prior_c, prior_n;
    double nprior_a, nprior_b;  //around 1
    double cprior_a, cprior_b;  //around ratio

    time_t last_time;
    unsigned long int seed;
    unsigned long int seed2;
    friend class GTree;


public:
    int N;
    int C;  //total number of elements
    int num_base;
    
    vector< string > species_names; // size S
    vector<string> nodes_names;  //size N
    // GSL random number generator
    gsl_rng * RNG;
    std::mt19937 twister;
    std::mt19937 twister2; //for shuffling input data only.
    
    // species tree (using array to accelerate)
    int    (*children)[2];
    int     *parent;
    double  *distances;
    double  *thetas;
    vector<double> heights;
    vector<int> move_br; //branches that may diff from Sp tree

    vector< vector< double >>  log_liks_WL;
    vector<vector< double >>   log_liks_Z; //Store max posterior 
    vector<vector<double>> log_mle;
    
    double br_sample_cutoff;
    
    //Han*: BPP dirichlet prior param arguments added
    //BPP(int pC, PhyloProf & _prof, PhyloTree & _tree, string output_path, string _target, string _outgroup, double _conserve_prop, string _conservegroup, double _ratio0, double _ratio1, int _ropt, double _cub, double _nlb, double _npriora, double _npriorb, double _cpriora, double _cpriorb, int _seed, double _prep_grate, double _prep_lrate, double _prep_lrate2, double _prior_g_a, double _prior_g_b, double _prior_l_a, double _prior_l_b,double _prior_l2_a, double _prior_l2_b, double _indel, double _indel2, double missing_thres, bool _sample_indel)
    BPP(int pC, PhyloProf & _prof, PhyloTree & _tree, string output_path, string _target, string _outgroup, double _conserve_prop, string _conservegroup, double _ratio0, double _ratio1, int _ropt, double _cub, double _nlb, double _npriora, double _npriorb, double _cpriora, double _cpriorb, int _seed, int _seed2, double _prep_grate, double _prep_lrate, double _prep_lrate2, double _prior_g_a, double _prior_g_b, double _prior_l_a, double _prior_l_b,double _prior_l2_a, double _prior_l2_b, double _indel, double _indel2, double missing_thres, bool _sample_indel, vector<double> _prior_dir_param, double _br_sample_cutoff, string _deepcoal_species)
{

    seed = _seed;
    seed2= _seed2;
    RNG = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(RNG, seed);
    
    twister.seed(seed);
    twister2.seed(seed2);

    C = _prof.element_names.size();
    G = _prof.G;
    S = _tree.S;
    N = S*2 - 1;

    //Han*:
    pi = _tree.pi;
    //for(int i=0; i<4; i++){
    //    pi[i]=_tree.pi[i];
    //}
    prior_dir_param=_prior_dir_param;

    ratio0 = _ratio0;
    ratio1 = _ratio1;
    cub = _cub;  // upper bound for conserved rate
    nlb = _nlb;  // lower bound for neutral rate
    ropt = _ropt; // restrict rate option: 0,1,2

    indel = _indel; //0.025; 0.25 for mammal
    indel2 = _indel2; //0.013 for ratite; 0.08 for mammal

    ind_lrate = _prep_lrate;
    ind_lrate2 = _prep_lrate2;
    ind_grate = _prep_grate;

    prior_l_a = _prior_l_a, prior_l_b = _prior_l_b;
    prior_l2_a = _prior_l2_a, prior_l2_b = _prior_l2_b;
    prior_g_a = _prior_g_a, prior_g_b = _prior_g_b;

    nprior_a = _npriora, nprior_b = _npriorb;  //around 1
    cprior_a = _cpriora, cprior_b = _cpriorb;  //around ratio

    // matching the species in profile and tree data
    MatchProfAndTree(_prof, _tree);

    element_size = vector<unsigned int>(C);
    element_start = vector<unsigned int>(C);
    element_tree = vector<string>(C);

    // initialization of profile and tree
    species_names = _prof.species_names;
    nodes_names = _tree.nodes_names; //first S elements matches with species_names.

    br_sample_cutoff = _br_sample_cutoff;

    if(!(_deepcoal_species=="")){ //user input constrained species.
        vector<string> tmp = strutils::split(strutils::trim(_deepcoal_species), ';');
        for (vector<string>::iterator it = tmp.begin(); it < tmp.end(); it++){
            ptrdiff_t pos = find(nodes_names.begin(), nodes_names.end(), *it) - nodes_names.begin();
            dp_species.push_back(pos);
        }
    }

    vector<string> tmp = strutils::split(strutils::trim(_target),';');
    for(vector<string>::iterator it = tmp.begin(); it<tmp.end();it++)
    {
        ptrdiff_t pos = find(species_names.begin(), species_names.end(), *it) - species_names.begin();
        target_species.push_back(pos);
    }

    tmp = strutils::split(strutils::trim(_conservegroup),';');
    for(vector<string>::iterator it = tmp.begin(); it<tmp.end();it++)
    {
        ptrdiff_t pos = find(species_names.begin(), species_names.end(), *it) - species_names.begin();
        conservedgroup.push_back(pos);
    }
    conserve_prop = _conserve_prop;

    if(indel<1e-10)
    {
        num_base = 4;
    }else{
        num_base = 5;

    }

    int start = 0;

    for(int c =0; c<C;c++)  // in mammal bed, the element pos is not continous
    {
        int end = start + _prof.element_pos[c][1] - _prof.element_pos[c][0];
        int len = end -start;
        element_size[c] = len;
        element_start[c] = start;
        //if(_prof.element_tree.size() >0) element_tree[c] =_prof.element_tree[c];
        if(_prof.element_tree.size() >0){
            element_tree[c] =_prof.element_tree[c]; //if input genetree, use input, and fixtree
        }
        start = end;
    }

    InitPhyloTree(_tree);  //, _indel here read in Q, depends on _indel!

    //reecord sp tree, in bppc.initMCMC pass on to initiate GTree
    
    if(_prof.element_tree.size() ==0){ 
        std::stringstream buffer;
        int rtN=N;
        for(int i=(N-1); i>=0; i--){
            if(parent[i]==N){
                rtN=i;
                break;
            }
        }

        getTreeString(rtN,buffer);
        element_tree[0]=buffer.str();
        for(int c=1; c<C; c++) element_tree[c]=element_tree[0];
    }


    cout << "InitPhyloTree finished" <<endl;

    // transition probability of hidden state Z
    TM_Int = vector<mat >(N, zeros<mat>(3,3));
    log_TM_Int = vector<mat >(N, zeros<mat>(3,3));

    tmp = strutils::split(strutils::trim(_outgroup),';');
    for(vector<string>::iterator it = tmp.begin(); it<tmp.end();it++)
    {
        ptrdiff_t pos = find(species_names.begin(), species_names.end(), *it) - species_names.begin();
        outgroup.push_back(pos);
    }
    getUppertree(N-1, outgroup, upper); //upper: set containing all nodes from roots to any species in outgroup
    getUppertree(N-1, conservedgroup, upper_conserve); //include root!

    if(num_base > 4)  // no indel
    {
        submat *= (1-indel);
        mat B = ones<mat>(4,1) * indel;
        mat C = ones<mat>(1,5) * indel2;

        submat.insert_cols(4, B);
        submat.insert_rows(4, C);

        colvec c = sum(submat,1);
        submat.diag() -=c;
        mat a = null(submat.t());
        pi = a/accu(a);
    }

    cx_mat bvec;
    cx_vec aval;
    eig_gen(aval, bvec, submat);
    eigenval = conv_to<mat>::from(aval);
    eigenvec = conv_to<mat>::from(bvec).t();
    eigeninv = inv(eigenvec);

    log_pi = log(pi);

    //Han*: calculate instantaneous rate
    inst_rate=zeros(6);
    inst_rate(0)=submat(0,1)/pi(1);
    inst_rate(1)=submat(0,2)/pi(2);
    inst_rate(2)=submat(0,3)/pi(3);
    //inst_rate(3)=submat(1,0)/pi(0);
    inst_rate(3)=inst_rate(0); //GTR;
    inst_rate(4)=submat(1,2)/pi(2);
    //inst_rate(5)=submat(1,3)/pi(3);
    inst_rate(5)=inst_rate[1];

    move_br = vector<int> (N, 0);
    cons_sample = false;
    if(dp_species.size()>0){ 
        for(int i=0; i<dp_species.size(); i++){
            bool upper_rel = true;
            set<int>:: iterator it_out=upper.find(dp_species[i]);
            if(it_out ==upper.end()){
                it_out = upper.find(parent[dp_species[i]]);
                if(it_out == upper.end()) upper_rel = false;
            }
            if((children[dp_species[i]][0] != -1) && (!upper_rel)){
                move_br[dp_species[i]] = 1;
                if(!cons_sample) cons_sample = true;
            } 
        }
        if(br_sample_cutoff != 10.0){
            for (int i = S; i < (N - 1); i++){
                bool upper_rel = true;
                set<int>:: iterator it_out=upper.find(i);
                if(it_out ==upper.end()){
                    it_out = upper.find(parent[i]);
                    if(it_out == upper.end()) upper_rel = false;
                }
                if (distances[i] < br_sample_cutoff && children[i][0] != -1 && (!upper_rel)){
                    move_br[i] = 1;
                    if(!cons_sample) cons_sample = true;
                }  
            }
        }
    }else{
        for (int i = S; i < (N - 1); i++){
            bool upper_rel = true;
            set<int>:: iterator it_out=upper.find(i);
            if(it_out ==upper.end()){
                it_out = upper.find(parent[i]);
                if(it_out == upper.end()) upper_rel = false;
            }
            if (distances[i] < br_sample_cutoff && children[i][0] != -1 && (!upper_rel)){
                //move_br[children[i][0]] = 1; //actual br whose gnode can be moved.
                //move_br[children[i][1]] = 1;
                move_br[i] = 1; //short br. later when init gene tree, based on this, put its child-gen to var_br_node.
                if((br_sample_cutoff<10.0) && (!cons_sample)) cons_sample = true;
            }
        }
    }
    //find the species node that is the root of all ingroup species. Top at this node and its sibling (i.e., along out group path), and above will not be touched.
    root_ingrp = N-1; //root
    for(int i = S; i<(N-1); i++){
        set<int>:: iterator it_out=upper.find(i);
        if(it_out==upper.end()){
            it_out=upper.find(parent[i]);
            if(it_out!=upper.end()){
                root_ingrp = i;
                break;
            }
        }
    }
    //cout<<"root node leading to all ingroup species is "<<root_ingrp<<endl;
    //for(set<int>:: iterator t1=upper.begin(); t1!=upper.end(); t1++){
    //    cout<<"outgrp "<< *t1 <<", sp="<< nodes_names[*t1]<<endl;
    //}
}

~BPP()
{

    delete [] children;
    delete [] parent;
    delete [] distances;
    delete [] thetas;

    gsl_rng_free(RNG);
}

static vec log_exp_multi(mat x, vec log_y)
{
    double log_ymax = arma::max(log_y);
    log_y = log_y - log_ymax;
        vec result = x*exp(log_y);

    return(log(result) + log_ymax);
}


static vec log_multi(mat log_x, vec log_y)  //X^T * y
{
    log_x.each_col()+=log_y;

    double log_xmax = log_x.max();
    if(log_xmax == -INFINITY)
    {
        vec xx(log_x.n_cols);
        xx.fill(-INFINITY);
        return(xx);
    }


    log_x -= log_xmax;

    rowvec result = log(sum(exp(log_x))) + log_xmax;

    return(result.t());
}

static vec log_multi2(mat log_x, vec log_y)
{
    log_x.each_col()+=log_y;
    for(size_t i =0; i < log_x.n_cols; i++)
    {
        double log_xmax = log_x.col(i).max();
        log_x.col(i) -= log_xmax;
        log_y.row(i) = log(sum(exp(log_x.col(i)))) + log_xmax;
    }
    return(log_y);
}

static rowvec log_exp_colsum(mat log_y)
{
    double log_ymax = log_y.max();
    log_y = log_y - log_ymax;

    rowvec result = log(sum(exp(log_y)));
    return(result + log_ymax);
}

static double log_exp_sum(mat log_y)
{
    double log_ymax = log_y.max();
    if(log_ymax == -INFINITY) return(-INFINITY);
    log_y = log_y - log_ymax;

    double result = accu(exp(log_y));
    return(log(result) + log_ymax);
}

static vec log_sample(vec log_y){  //return unnormalized prob
    double log_ymax = arma::max(log_y);
    log_y -= log_ymax;
    return(exp(log_y));
}

static vec log_sample_norm(mat log_y){  //return normalized prob
        double x = log_exp_sum(log_y);
    if(x == INFINITY){ log_y.fill(0); return(log_y);}
        log_y -= x;
        log_y = exp(log_y);
    return(log_y);
}

static void printZ(int N, vector<int> & Z, int (*children)[2])
{
    bool precol =false;
    queue<int> stack;
    vector<bool> colors = vector<bool> (N,false);

    stack.push(N-1);

    while(stack.size()>0)
    {
        int j = stack.front();

        if(colors[j]!=precol) cout<<endl;
        cout << j << ": "<< Z[j] <<", ";

        precol = colors[j];


        stack.pop();

        if(children[j][0]!=-1)
        {
            for(int chi=0;chi<2;chi++)
            {
                int c = children[j][chi];
                colors[c]=!colors[j];
                stack.push(c);

            }

        }
    }
    cout<<endl;
    cout<<endl;

}

void MatchProfAndTree(PhyloProf & _prof, PhyloTree & _tree);
void InitPhyloTree(PhyloTree & tree);//double indel, double indel2 = 0.01, double indel_pi = 0.5);
void InitMCMC(int _num_burn, int _num_mcmc, int _num_thin);
void Output_init(PhyloProf & prof, string output_path, vector<int> & ids);
void Output_init0(PhyloProf & prof, ofstream& out_lik, vector<int> & ids);
void Output_simu(PhyloProf & prof, string outpath, int i);
   
void getSubtree(int root, vector<int> & visited_init);
void getSubtree(int root, set<int>& child, vector<int> & visited_init);
void getUppertree(int root, vector<int>& child, set<int> & visited_init);

void sample_indel(vector< string> & X, double missing_thres, string output_path, int _iter = 100, int _MH = 10, bool _verbose = true);
void sample_proposal(int iter, double & lrate_prop, double & grate_prop, ofstream & output);
void sample_hyperparam(int iter, vector<int> & ids, ofstream & output); //double indel_prop, double indel2_prop,
vec getlogTM(double dist, double rate, mat & lam );
mat getlogTM(double dist, double rate);
//Han*: declare getlogTMc function
mat getlogTMc(double dist, double rate, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv);
double log_lik(vector< vector<vec> > & lambda, double indel, double indel2, int iter, int block, vector<unsigned int> & v, double p);

mat getlogTM_len(double dist, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv); //25-May for Sample_tree2, directly input total r*t.

void getTreeString(int rootN, std::stringstream & buffer);
};

#endif /* bpp_hpp */
