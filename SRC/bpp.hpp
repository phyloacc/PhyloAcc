//
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

//const double T1 = 0.5;

const double LOG_ZERO = -INFINITY;

class BPP
{
    friend class BPP_C;
    
private:
    
    int G;  //total base pairs
    int S;
    vector <unsigned int> element_size;
    vector <unsigned int> element_start;
    
    double ratio0;
    double ratio1;
    double cub;
    double nlb;
    int ropt;
    
    // names of species
    vector<int> target_species;
    vector<int> conservedgroup;
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
    
    
    // phylogenetic tree (using array to accelerate)
    int    (*children)[2];
    int     *parent;
    double  *distances;
    vec pi;
    vec log_pi;
    mat eigenvec; vec eigenval;
    mat eigeninv;
    mat submat;
    double indel;
    double indel2;
    
    int num_base;
    
    mat eigenvecprop; vec eigenvalprop;
    mat eigeninvprop; vec piprop;
    
    int moveroot = -10; // move root to one of its children
    
    
    vector<vector < vector <int > >> Max_Z;
    vector<vector < vector <int > >> cur_Z;
    
    
    double ind_lrate, ind_lrate2;  //for Z
    double ind_grate;
    double vlr = 100; // proposal variance of loss rate
    double vgr = 100; // proposal variance of gain rate
    //double prior_glr[3];
    double prior_l_a, prior_l_b;
    double prior_l2_a, prior_l2_b;
    double prior_g_a, prior_g_b;
    
    vector<mat >  TM_Int;
    //vector<mat>  TM_null;
    vector<mat> log_TM_Int;
    
    
    //vector<mat > log_cache_TM_neut;
    //vector<mat > log_cache_TM_cons;
    vector<mat > log_cache_TM_null;
    
    vector< double >  log_liks_sgl;
    vector<vector< double >>   log_liks_Z; //Store max logliklihood
    vector< double > log_liks_curZ; //Store current loglik
    vector< double > log_liks_propZ; //Store loglik with prop indel parameter
    
    vector< double > MH_ratio_gain; //
    vector< double > MH_ratio_loss; //
    
    vector< double > log_liks_null;
    vector< double > log_liks_resZ;
    
    vector< double > cur_nrate;
    vector< double > cur_crate;
    
    vector< double > cur_lrate;
    vector< double > cur_grate;
    vector< double > cur_lrate2;

    
    vector < vector <mat>> log_cache_TM_neut_ensemble;
    vector < vector <mat>> log_cache_TM_cons_ensemble;
    vector<double> neut_prior_density;
    vector<double> cons_prior_density;
    
    
    //hyperparameter of prior_c, prior_n;
    double nprior_a, nprior_b;  //around 1
    double cprior_a, cprior_b;  //around ratio
    
    
    
    
    time_t last_time;
    unsigned long int seed;


    
public:
    int N;
    int C;  //total number of elements
    vector< string > species_names; // size S
    vector<string> nodes_names;  //size N
    // GSL random number generator
    gsl_rng * RNG;



    BPP(int pC, PhyloProf & _prof, PhyloTree & _tree, string output_path, string _target, string _outgroup, double _conserve_prop, string _conservegroup, double _ratio0, double _ratio1, int _ropt, double _cub, double _nlb, double _npriora, double _npriorb, double _cpriora, double _cpriorb, int _seed, double _prep_grate, double _prep_lrate, double _prep_lrate2, double _prior_g_a, double _prior_g_b, double _prior_l_a, double _prior_l_b,double _prior_l2_a, double _prior_l2_b, double _indel, double _indel2, double missing_thres, bool _sample_indel)
{
    
    seed = _seed;
    RNG = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(RNG, seed);

    
    C = _prof.element_names.size();
    G = _prof.G;
    S = _tree.S;
    N = S*2 - 1;
    pi = _tree.pi;
    
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
//    prior_glr[0] = _prior_l * _prior_gla;
//    prior_glr[1] = _prior_g * _prior_gla;
//    prior_glr[2] = (1 - _prior_l -_prior_g)  * _prior_gla;
    
    prior_l_a = _prior_l_a, prior_l_b = _prior_l_b;
    prior_l2_a = _prior_l2_a, prior_l2_b = _prior_l2_b;
    prior_g_a = _prior_g_a, prior_g_b = _prior_g_b;
    
    
    
    nprior_a = _npriora, nprior_b = _npriorb;  //around 1
    cprior_a = _cpriora, cprior_b = _cpriorb;  //around ratio
    
    //nprior_a = 15, nprior_b = 0.1;  //around 1
    //cprior_a = 5, cprior_b = 0.02;  //around ratio
    
    
    // matching the species in profile and tree data
    MatchProfAndTree(_prof, _tree);
    
    element_size = vector<unsigned int>(C);
    element_start = vector<unsigned int>(C);
    
    // initialization of profile and tree
    species_names = _prof.species_names;
    nodes_names = _tree.nodes_names;
    
    //_refspecies = strutils::trim(_refspecies);
    //refspecies = find(species_names.begin(), species_names.end(),  _refspecies) - species_names.begin();
    
    
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
    
    //X = vector < vector <string> >(C);
    
    int start = 0;
    
    for(int c =0; c<C;c++)  // in mammal bed, the element pos is not continous
    {
        //int start = _prof.element_pos[c][0];
        //int end = _prof.element_pos[c][1];
        int end = start + _prof.element_pos[c][1] - _prof.element_pos[c][0];
        int len = end -start;
        element_size[c] = len;
        element_start[c] = start;
//        for(int s=0; s<S;s++){
//            //string y =strutils::ToLowerCase(_prof.X[s].substr(start,end));  // already to lower in reading alignment
//            X[c].push_back(_prof.X[s].substr(start,end));
//        }
        
        start = end;
        
    }
    

    InitPhyloTree(_tree);  //, _indel here read in Q, depends on _indel!
    
    //cout <<"eigenval: " << eigenval;
    //cout <<"eigenvec: " << eigenvec;
    
    cout << "InitPhyloTree finished" <<endl;
    
    // transition probability of hidden state Z
    TM_Int = vector<mat >(N, zeros<mat>(3,2));
    log_TM_Int = vector<mat >(N, zeros<mat>(3,2));
    
    
    
//    for(int s=0; s<N; s++)
//    {
//
//        //double x = exp(-ind_grate *distances[s]);
//        TM_Int[s](0,0) = 1 - ind_grate;
//        TM_Int[s](1,0) = ind_grate * (1 - ind_lrate2);
//        TM_Int[s](2,0) = ind_grate * ind_lrate2;
//        //double y = exp(-ind_lrate *distances[s]);
//        TM_Int[s](1,1) = 1 - ind_lrate;
//        TM_Int[s](2,1) = ind_lrate;
//
//    }

    tmp = strutils::split(strutils::trim(_outgroup),';');
    for(vector<string>::iterator it = tmp.begin(); it<tmp.end();it++)
    {
        ptrdiff_t pos = find(species_names.begin(), species_names.end(), *it) - species_names.begin();
        outgroup.push_back(pos);
    }
    getUppertree(N-1, outgroup, upper);
    getUppertree(N-1, conservedgroup, upper_conserve);
    
    for(int s=0; s<N-1; s++) // nodes: from bottom to top, doesn't include root
    {
        if(upper.find(s)==upper.end())
        {
            subtree.push_back(s);
        }
    }

    
   // ctnutils::DispVector(upper);
   // cout << endl;
    
//    for(set<int>:: iterator it = upper.begin(); it!=upper.end();it++)
//    {
//        TM_Int[*it](1,1) = 1;
//        TM_Int[*it](2,1) = 0;
//
//        TM_Int[*it](0,0) = 1 - ind_grate;
//        TM_Int[*it](1,0) = ind_grate;
//        TM_Int[*it](2,0) = 0;
//    }
//
//
//    for(int s=0; s<N; s++){
//        log_TM_Int[s] = log(TM_Int[s]);
//    }
    
    
    //if(_sample_indel) sample_indel(_prof.X, missing_thres, output_path, 200, 25); //_prof.X
    
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
    
    log_cache_TM_null = vector<mat >(N, zeros<mat> (num_base,num_base));
    
    
    
    for(int i=0;i<N-1;i++)
    {
        
        if(distances[i]> 1e-10 )
        {
            mat tmp_diag  = exp(eigenval*distances[i]);
            mat x = eigenvec;
            x.each_col()%=tmp_diag;
            
            //log_cache_TM_neut[i] = log(eigeninv * x); //transpose Q
            log_cache_TM_null[i] = log(eigeninv * x); //log_cache_TM_neut[i] ;
            
//            tmp_diag  = exp(eigenval*distances[i]*ratio); //*scale_cons[c]);
//            x = eigenvec;
//            x.each_col()%=tmp_diag;
//            
//            log_cache_TM_cons[i] = log(eigeninv * x);
            
        }else{
            
            log_cache_TM_null[i].fill(-INFINITY); //83, root
            log_cache_TM_null[i].diag().fill(0);
            //log_cache_TM_cons[i] = log_cache_TM_neut[i] ;
            //log_cache_TM_null[i] = log_cache_TM_null[i] ;
            //log_cache_TM_neut[i] = log_cache_TM_neut[i] ;
        }
        
    }
    
    
        
    //cout << "cache finished" <<endl;

   
}

~BPP()
{
    
    delete [] children;
    delete [] parent;
    delete [] distances;
    
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
    
   for(size_t i =0; i < log_x.n_cols; i++){
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
    
static vec log_sample_norm(mat & log_y){  //return normalized prob
        double x = log_exp_sum(log_y);
    if(x == INFINITY){ log_y.fill(0); return(log_y);}
        log_y -= x;
        log_y = exp(log_y);
    return(log_y);
}
    
static void printZ(int N, vector<int> & Z, int (*children)[2])
{
    //int j = root;
    
    //cout << nodes_names[j]<<"\t";
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
void getSubtree(int root, vector<int> & visited_init);
void getSubtree(int root, set<int>& child, vector<int> & visited_init);
void getUppertree(int root, vector<int>& child, set<int> & visited_init);

    void sample_indel(vector< string> & X, double missing_thres, string output_path, int _iter = 100, int _MH = 10, bool _verbose = true);
void sample_proposal(int iter, double & lrate_prop, double & grate_prop, ofstream & output);
void sample_hyperparam(int iter, vector<int> & ids, ofstream & output); //double indel_prop, double indel2_prop, 
double log_lik(vector< vector<vec> > & lambda, double indel, double indel2, int iter, int block, vector<unsigned int> & v, double p);
    };

#endif /* bpp_hpp */
