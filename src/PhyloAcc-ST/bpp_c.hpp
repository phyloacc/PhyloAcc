//
//  bpp_c.hpp for mammal
//  PhyloAcc_init2
//
//  Created by hzr on 4/19/16.
//  Copyright © 2016 hzr. All rights reserved.
//

#ifndef bpp_c_hpp
#define bpp_c_hpp


#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <string>
#include <armadillo>
#include <cassert>

#include <cmath>
#include <sys/time.h>

#include <map>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#include <memory>

#include "../PhyloAcc-common/newick.h"
#include "../PhyloAcc-common/profile.h"
#include "../PhyloAcc-common/bpp_constructor.h"
#include "../PhyloAcc-common/rng.h"
#include "bpp.hpp"
#include "../PhyloAcc-common/utils.h"


using namespace std;
using namespace arma;

class BPP_C
{
    private:
    
    int CC = 0;  //current number of elements
    int GG = 0; //base pairs current elements
    
    int S = 0;
    int N = 0;
    
    vector< vector<vec > > ambiguousS_null;  //base * species# * 4;
    
    std::unique_ptr<int[][2]> children2;  //children from pruned tree
    std::unique_ptr<double[]> distances2;
    std::unique_ptr<int[]> parent2;
    int root = -1;
    vec pi;
    vec log_pi;
    
    
    vector<int> nodes;
    vector<int> internal_nodes;
    vector<bool> missing;
    vector<int> upper_c;
    vector<int> upper_conserve_c;
    
   // double prior_glr[3];
    double prior_l_a = 0.0, prior_l_b = 0.0;
     double prior_l2_a = 0.0, prior_l2_b = 0.0;
     double prior_g_a = 0.0, prior_g_b = 0.0;
    double ratio0 = 0.0;
    double ratio1 = 0.0;
    int num_burn = 0;   // num of burn-in updates
    int num_mcmc = 0;   // num of MCMC updates
    int num_thin = 0;   // num of updates between two samples
    int adaptive_freq = 100;
   
    vector<int>Z ; //N *0, accelerate(time,0-1), loss(-1)
    vector<int>fixZ ;
    vector < vector<int> > Tg; //current history for each element, GG*N
    
    
    
    //buffer P(zt|zt-1,Hg)
    vector<vec> log_prob_back;
    
    //get P(X|theta)
    vector< vector<double> >  log_emission; //( N-1)*2
    
    
    
    // MCMC updating states
    int m = 0;                                      // current MCMC step
    
    vector<mat> log_TM_Int;
    
    vec prior_z;
    
    vector<mat > log_cache_TM_neut;
    vector<mat > log_cache_TM_cons;
    vector<mat > log_cache_TM_null;
    
    
    vector< vector<vec> > lambda;  //GG *N *acgt
    

    
    
    // samples to output
    double MaxLoglik = -INFINITY;
    int Max_m = 0;
    vector <int > Max_Z;
    
    vector< double >  trace_loglik;  //P(X|Z, TM, r)
    vector< double >  trace_full_loglik;  //P(X,Z,r|TM)
    vector< vector<int> > trace_Z;  //iteration * Species * {0,1,2}, only for init
    vector <double> trace_n_rate;  //iter*element
    vector <double> trace_c_rate;
    
    vector <double> trace_l_rate;  //iter*element
    vector <double> trace_l2_rate;
    vector <double> trace_g_rate;
    
    int accept_n_rate = 0;
    int accept_c_rate = 0;  //how many accepted in current cycle
    
    double prop_n = 0.0;  //for adaptive MCMC, changed by acceptance rate
    double prop_c = 0.0;
    
    
    
    // GSL random number generator
    gsl_rng * RNG = nullptr;
    
    time_t last_time = 0;
    
    unsigned long int seed = 0;

    
public:
    BPP_C(const BPP_C&) = delete;
    BPP_C& operator=(const BPP_C&) = delete;

    bool failure = false;
    bool verbose = false; 
    
    BPP_C(int c, PhyloProf _prof, BPP& bpp, char gapchar, double missing_thres, bool & filter, bool _verbose, double consToMis, bool prune=0, double revgap=0, int min_length =50, double nconsToMis = 1, int chain_index = 0)//, double _indel)
    {
        
        RNG = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(RNG, phyloacc::DeriveSeed(bpp.seed, phyloacc::ProgramKind::ST,
                                              chain_index, c, 0,
                                              phyloacc::RngStream::WorkerGsl));
        
        num_burn = bpp.num_burn;   // num of burn-in updates
        num_mcmc = bpp.num_mcmc;   // num of MCMC updates
        adaptive_freq = bpp.num_thin;   // num of updates between two samples, 100
        
       
        N = bpp.N;
        CC = c;
        GG =bpp.element_size[c];

        if(GG < min_length)
        {
            filter = true;
            return;
        }
        
        root = N-1;
        S = bpp.S;
        
        
        ratio0 = bpp.ratio0;  //no use
        ratio1 = bpp.ratio1; //no use
        
        verbose = _verbose;

        
        log_cache_TM_neut = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_cons = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_null = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        
        log_TM_Int = vector<mat >(N, zeros<mat>(3,3));
        
//        for(size_t i =0 ;i <3; i++)  // hyperparameter for loss and gain rates
//        {
//            prior_glr[i] = bpp.prior_glr[i];
//        }
       
        
        prior_l_a = bpp.prior_l_a;
        prior_l_b = bpp.prior_l_b;
        prior_g_a = bpp.prior_g_a;
        prior_g_b = bpp.prior_g_b;
        prior_l2_a = bpp.prior_l2_a;
        prior_l2_b = bpp.prior_l2_b;
        
        
        prior_z = zeros<vec>(3);
        prior_z[0] = 0.5;  // prior for root
        prior_z[1] = 0.5;
        prior_z = log(prior_z);
        
        
        if(verbose) cout << "Init lambda" <<endl;
        
        int st = bpp.element_start[CC];
        vector<int> site_order;
        phyloacc::LeafEncoding leaf_encoding = phyloacc::EncodeLeafAlignment(
            _prof.X, st, GG, S, N, bpp.num_base, gapchar,
            phyloacc::MissingBasePolicy::GapOnly, site_order);
        lambda = leaf_encoding.lambda;
        Tg = leaf_encoding.tg;

        phyloacc::ColumnFilterResult column_filter = phyloacc::RemoveHighMissingColumns(
            lambda, Tg, S, revgap, min_length);
        if(column_filter.filtered)
        {
            filter = true;
            return;
        }
        GG = column_filter.length;

        vector<int> num_missing = phyloacc::CountMissingBySpecies(Tg, S);

        ambiguousS_null = vector<vector<vec> > (GG, vector<vec>(S,zeros<vec>(bpp.num_base)));
        
        
        children2.reset(new int[N][2]);
        parent2.reset(new int[N]);
        distances2.reset(new double[N]);
        for(int i=0; i<N; i++)
        {
            children2[i][0] = bpp.children[i][0];
            children2[i][1] = bpp.children[i][1];
            parent2[i] = bpp.parent[i];
            distances2[i] = bpp.distances[i];
        }
       
        // set missing
        missing = phyloacc::BuildMissingNodes(S, N, children2.get(), num_missing, missing_thres, GG);
        
//        cout << "missing: ";
//        for(int s = 0; s<N;s++)
//        {
//            if(missing[s]) cout << s << " ";
//        }
//        cout <<endl;
        
        
        if(phyloacc::ConservedMissingExceeds(missing, bpp.conservedgroup, bpp.conserve_prop))
        {
            filter = true;
            return;
        }
        
        
        
        
        //prune tree if outgroup not conserved, find root
        if(prune)
        {
            set<int> alls;
            for(int s=0; s<N; s++) alls.insert(s);
            getSubtree_missing(N-1, alls, -1);
        }else{
            getSubtree_missing(N-1, bpp.upper, -1);
        }
        
        
        
        parent2[root] = N;  //parent2 not correct for all nodes, but children2 is correct
    
        if(verbose) cout << "root: " << root << endl;
        
        
        getSubtree(root, nodes);
        for(vector<int>::iterator it = nodes.begin(); it < nodes.end(); it++)
        {
            if(*it>=S)
            {
                internal_nodes.push_back(*it);  // internal nodes
            }

        }
        phyloacc::CollectUpperNodesInSubtree(nodes, bpp.upper, bpp.upper_conserve,
                                             upper_c, upper_conserve_c);
        
        
        // if both children are missing, can't infer the parent's base pair, just set to 'missing'
        for(vector<int>::iterator it = internal_nodes.begin(); it < internal_nodes.end(); it++)
        {
            int* p = children2[*it];
            for(int g=0; g<GG; g++){
                if(Tg[g][p[0]] >= bpp.num_base && Tg[g][p[1]]>= bpp.num_base)
                {
                    Tg[g][*it] = bpp.num_base;
                    
                }
            }
            
        }
        
        log_pi = bpp.log_pi;
        
  
        // from initMCMC
        phyloacc::BppCTraceBuffers traces = phyloacc::InitializeBppCTraceBuffers(
            num_burn+num_mcmc, N, 0, 0, 0);
        trace_loglik = traces.trace_loglik; //P(X|Z, r)
        trace_full_loglik = traces.trace_full_loglik; //P(X, Z, r)
        trace_Z = traces.trace_z;
        trace_n_rate = traces.trace_n_rate;
        trace_c_rate = traces.trace_c_rate;
        trace_l_rate = traces.trace_l_rate;
        trace_g_rate = traces.trace_g_rate;
        trace_l2_rate = traces.trace_l2_rate;

        //getEmission_ambig();
        
        
        
        
        log_emission = traces.log_emission;
        
        
        for(int s=0;s<N;s++)  // For all nodes !!! ....only terminal nodes, S
        {
            if(missing[s])
            {
                //            log_prob_back[s][2] += log(nconsTomis);
                //            log_prob_back[s][1] = log(consTomis) ;//-INFINITY;
                //            log_prob_back[s][0] += log(nconsTomis);
                
                log_emission[s][2] = log(nconsToMis);
                log_emission[s][1] = log(consToMis) ;//-INFINITY;
                log_emission[s][0] = log(nconsToMis);
                
                //fixZ[s] = 0; // ?
                
            }
            //        }else{
            //            
            //            log_prob_back[s][2] += log(1 - nconsTomis);
            //            log_prob_back[s][1] = log( 1 - consTomis) ;//-INFINITY;
            //            log_prob_back[s][0] += log(1 - nconsTomis);
            //            
            //
            //        }
        }

    }
    
    
    
    
    ~BPP_C()
    {
        gsl_rng_free(RNG);
    }

    void getSubtree(int root, set<int>& child, vector<int> & visited_init);
    void getSubtree(int root, vector<int> & visited_init);
    void getSubtree_missing(int root, set<int>& upper, int child);
    
    void getEmission_ambig();
    void initMCMC(int iter, BPP&bpp, int resZ);
    void Update_Tg(int g, vector<bool> visited,BPP& bpp, bool tosample);
    void getUpdateNode(vector<int> changedZ, vector<bool> & visited_init);
    void MonitorChain(int m, BPP &bpp, double &loglik, const double add_loglik, const int resZ) ;
    void getUpdateNode(bool neut,vector<bool> & visited_init);
    double log_f_Xz(vector<bool> visited, vector<vector<vec>>& lambda_tmp, bool neut, double propos, BPP& bpp);
    double sample_rate(int resZ, double old_rate, bool neut, vector<bool> visited, double & loglik_old, BPP& bpp, int M =1, bool adaptive = true, double adaptive_factor = 0.5);
    void Gibbs(int iter, BPP &bpp, ofstream & outZ, string output_path,string output_path2,int resZ, bool UpR, bool UpHyper, double lrate_prop, double grate_prop);
    vector<int> Update_Z_subtree(int num_base = 5, bool prior = false);
    string Output_init_row(BPP& bpp, int resZ);
    void Output_init(string output_path, string output_path2, BPP& bpp,ofstream& outZ, int resZ);
    void Output_sampling(int iter, string output_path2, BPP &bpp, int resZ);
    
    void getEmission(int num_base);
    
    vector<int>  Move_Z(int & propConf, int & revConf,  int & changeZ);
    double log_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut,  vector<mat> & log_cache_TM_cons);
    double Update_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut, vector<mat> & log_cache_TM_cons, vector<bool> & visited);
    void log_f_Z(vector<int>& Z, vector<mat> & log_Int, double & MH_ratio_g, double & MH_ratio_l);
    
    //void Eval(BPP&bpp,int resZ); //, int numH,int numHZ);
    void Eval2(BPP&bpp, int resZ);
    //double prior_Z_subtree(vector<int> & tmpZ) ;
    vector<double> prior_Z_subtree(vector< vector<int> > & configZ, vector< int > numConfigZ);
    void sample_transition( double  & gr, double  & lr, double  & lr2);
    
};


#endif /* bpp_c_hpp */
