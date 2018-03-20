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

#include "newick.h"
#include "profile.h"
#include "bpp.hpp"
#include "utils.h"


using namespace std;
using namespace arma;

class BPP_C
{
    private:
    
    int CC;  //current number of elements
    int GG; //base pairs current elements
    
    int S;
    int N;
    
    vector< vector<vec > > ambiguousS_null;  //base * species# * 4;
    
    int    (*children2)[2];  //children from pruned tree
    double  *distances2;
    int *parent2;
    int root;
    vec pi;
    vec log_pi;
    
    
    vector<int> nodes;
    vector<int> internal_nodes;
    vector<bool> missing;
    vector<int> upper_c;
    
    double ratio0;
    double ratio1;
    int num_burn;   // num of burn-in updates
    int num_mcmc;   // num of MCMC updates
    int num_thin;   // num of updates between two samples
    int adaptive_freq = 100;
   
    vector<int>Z ; //N *0, accelerate(time,0-1), loss(-1)
    vector<int>fixZ ;
    vector < vector<int> > Tg; //current history for each element, GG*N
    
    
    
    //buffer P(zt|zt-1,Hg)
    vector<vec> log_prob_back;
    
    //get P(X|theta)
    vector< vector<double> >  log_emission; //( N-1)*2
    
    
    
    // MCMC updating states
    int m;                                      // current MCMC step
    
    vector<mat> log_TM_Int;
    
    
    vector<mat > log_cache_TM_neut;
    vector<mat > log_cache_TM_cons;
    vector<mat > log_cache_TM_null;
    
    
    vector< vector<vec> > lambda;  //GG *N *acgt
    

    
    
    // samples to output
    double MaxLoglik;
    int Max_m;
    vector <int > Max_Z;
    
    vector< double >  trace_loglik;  //P(X|Z, TM, r)
    vector< double >  trace_full_loglik;  //P(X,Z,r|TM)
    vector< vector<int> > trace_Z;  //iteration * Species * {0,1,2}, only for init
    vector <double> trace_n_rate;  //iter*element
    vector <double> trace_c_rate;
    
    int accept_n_rate = 0;
    int accept_c_rate = 0;  //how many accepted in current cycle
    
    double prop_n;  //for adaptive MCMC, changed by acceptance rate
    double prop_c;
    
    
    
    // GSL random number generator
    gsl_rng * RNG;
    
    time_t last_time;
    
    unsigned long int seed;
    bool verbose; //= false;

    
public:
    bool failure = false;
   
    
    BPP_C(int c, PhyloProf _prof, BPP& bpp, char gapchar, double missing_thres, bool & filter, bool _verbose, double consToMis, double nconsToMis = 1)//, double _indel)
    {
        
        RNG = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(RNG, bpp.seed);
        
        num_burn = bpp.num_burn;   // num of burn-in updates
        num_mcmc = bpp.num_mcmc;   // num of MCMC updates
        adaptive_freq = bpp.num_thin;   // num of updates between two samples, 100
        
       
        N = bpp.N;
        CC = c;
        GG =bpp.element_size[c];
        Tg = vector<vector<int> > (GG, vector<int>(N, -1));
        root = N-1;
        S = bpp.S;
        
        lambda = vector< vector<vec> >(GG, vector<vec>(N, zeros<vec>(bpp.num_base)));
        ambiguousS_null = vector<vector<vec> > (GG, vector<vec>(S,zeros<vec>(bpp.num_base)));
        
        
        ratio0 = bpp.ratio0;  //no use
        ratio1 = bpp.ratio1; //no use
        
        verbose = _verbose;

        
        log_cache_TM_neut = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_cons = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_null = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        
        log_TM_Int = vector<mat >(N, zeros<mat>(3,2)); 
        
        
        if(verbose) cout << "Init lambda" <<endl;
        vector<int> num_missing = vector<int> (S,0);
        
        int st = bpp.element_start[CC];
        //get leaves lambda
        for(int s=0; s<S; s++)
        {
            //const char* y = bpp.X[c][s].c_str();
            string y = _prof.X[s].substr(st, st+GG);
            
            for(int g=0; g<GG; g++){
                lambda[g][s].fill(LOG_ZERO);
                if( y[g]== gapchar) //y[g]=='n' ||
                {
                    num_missing[s] ++;
                    //Tg[g][s] = 4;  // for missing base pair
                }
                switch (y[g])
                {
                    case 'a':
                        lambda[g][s][0] = 0; Tg[g][s] = 0; break;
                    case 'c':
                        lambda[g][s][1] = 0; Tg[g][s] = 1; break;
                    case 'g':
                        lambda[g][s][2] = 0; Tg[g][s] = 2; break;
                    case 't':
                        lambda[g][s][3] = 0; Tg[g][s] = 3; break;
                    case 'r':
                        lambda[g][s][0] = 0;
                        lambda[g][s][2] = 0;
                        break;
                    case 'y':
                        lambda[g][s][1] = 0;
                        lambda[g][s][3] = 0;   break;
                    case 'k':
                        lambda[g][s][2] = 0;
                        lambda[g][s][3] = 0;  break;
                    case 'm':
                        lambda[g][s][0] = 0;
                        lambda[g][s][1] = 0;  break;
                    case 's':
                        lambda[g][s][1] = 0;
                        lambda[g][s][2] = 0;   break;
                    case 'w':
                        lambda[g][s][0] = 0;
                        lambda[g][s][3] = 0;   break;
                    case '-':
                        if(bpp.num_base <= 4){
                            for(int b =0;b<bpp.num_base;b++) lambda[g][s][b] = 0;
                        }
                        else{
                            lambda[g][s][4] = 0;
                        }
                        Tg[g][s] = 4; break;
                    default:
                        for(int b =0;b<bpp.num_base;b++) lambda[g][s][b] = 0;
                        Tg[g][s] = bpp.num_base;
                }
            }
            
        }

        
        
        children2    = new int[N][2];
        parent2      = new int[N];
        distances2   = new double[N];
        for(int i=0; i<N; i++)
        {
            children2[i][0] = bpp.children[i][0];
            children2[i][1] = bpp.children[i][1];
            parent2[i] = bpp.parent[i];
            distances2[i] = bpp.distances[i];
        }
       
        // set missing
        missing = vector<bool>(N, false);
        
        for(int s=S; s<N; s++)
        {
            int* p = children2[s];
            for(int cc=0;cc<2;cc++)
            {
                int chi = p[cc];
                if(chi<S && num_missing[chi]>missing_thres*GG)
                {
                    
                    missing[chi] = true;
                    
                }
            }
            
            if(missing[p[0]] && missing[p[1]]) {
                missing[s] = true;
                
            }
        }
        
//        cout << "missing: ";
//        for(int s = 0; s<N;s++)
//        {
//            if(missing[s]) cout << s << " ";
//        }
//        cout <<endl;
        
        
        int ct =0;
        for(vector<int>::iterator it = bpp.conservedgroup.begin(); it<bpp.conservedgroup.end();it++)
        {
            if(missing[*it]){
                ct +=1;
            }
        }
        
        if(ct > bpp.conserve_prop * bpp.conservedgroup.size())
        {
            filter = true;
            return;
        }
        
        
        
        
        //prune tree if outgroup not conserved, find root
        getSubtree_missing(N-1, bpp.upper, -1);
        parent2[root] = N;  //parent2 not correct for all nodes, but children2 is correct
    
        if(verbose) cout << "root: " << root << endl;
        

        
        getSubtree(root, nodes);
        for(vector<int>::iterator it = nodes.begin(); it < nodes.end(); it++)
        {
            if(*it>=S)
            {
                internal_nodes.push_back(*it);  // internal nodes
            }
            
            if(bpp.upper.find(*it)!=bpp.upper.end())
            {
                    upper_c.push_back(*it);
            }

        }
        
        
        // if both children are missing, can't infer the parent's base pair, just set to 'missing'
        for(vector<int>::iterator it = internal_nodes.begin(); it < internal_nodes.end(); it++)
        {
            int* p = children2[*it];
            for(int g=0; g<GG; g++){
                if(Tg[g][p[0]]== bpp.num_base && Tg[g][p[1]]==bpp.num_base)
                {
                    Tg[g][*it] = bpp.num_base;
                    
                }
            }
            
        }
        
        log_pi = bpp.log_pi;
        
  
        // from initMCMC
        trace_loglik = vector<double >(num_burn+num_mcmc, 0); //P(X|Z, r)
        trace_full_loglik = vector<double >(num_burn+num_mcmc, 0); //P(X, Z, r)
        trace_Z = vector<vector<int> >(num_burn+num_mcmc, vector<int>(N,1));
        trace_n_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_c_rate = vector<double >(num_burn+num_mcmc, 0);

        //getEmission_ambig();
        
        
        prop_c = 1.2; //(double) 1/GG;
        prop_n = 1.2; //(double) 1/GG * 10;
        
        
        log_emission = vector<vector<double> >(N-1, vector<double>(3,0));
        
        
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
    double sample_rate(int resZ, double old_rate, bool neut, vector<bool> visited, double & loglik_old, BPP& bpp, int M =1, bool adaptive = true, double adaptive_factor = 0.1);
    void Gibbs(int iter, BPP &bpp, ofstream & outZ, string output_path,string output_path2,int resZ, bool UpR, bool UpHyper, double lrate_prop, double grate_prop);
    vector<int> Update_Z_subtree(int num_base = 5);
    void Output_init(string output_path, string output_path2, BPP& bpp,ofstream& outZ, int resZ);
    void Output_sampling(int iter, string output_path2, BPP &bpp, int resZ);
    
    void getEmission(int num_base);
    
    vector<int>  Move_Z(int & propConf, int & revConf,  int & changeZ);
    double log_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut,  vector<mat> & log_cache_TM_cons);
    double Update_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut, vector<mat> & log_cache_TM_cons, vector<bool> & visited);
    void log_f_Z(vector<int>& Z, vector<mat> & log_Int, double & MH_ratio_g, double & MH_ratio_l);
    
    //void Eval(BPP&bpp,int resZ); //, int numH,int numHZ);
    void Eval2(BPP&bpp, int resZ);
    double prior_Z_subtree(vector<int> & tmpZ) ;
    
    
};


#endif /* bpp_c_hpp */
