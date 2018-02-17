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
    
    //vector< vector<vec > > ambiguousS_null;  //base * species# * 4;
    
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
    int adaptive_freq =  100;
   
    vector<int>Z ; //N *0, accelerate(time,0-1), loss(-1)
    vector<int>fixZ ;
    //vector < vector<int> > Tg; //current history for each element, GG*N
    //vector<double> B;
    vector<mat > posterior_Tg;
    vector<mat > posterior_Tg_next;
    
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
    
    vector<mat > log_cache_TM_neut_B;
    vector<mat > log_cache_TM_cons_B;
    vector<mat > log_cache_TM_null_B;
    
    vector< vector<vec> > lambda;  //GG *N *acgt
    vector< vector<vec> > message;  //GG *N *acgt
    
    //vector<vec> mean_lambda;  //N *acgt
    //vector<vec> mean_message;  //N *acgt
    
    vector<bool> isGB;
    
    // samples to output
    double logpZC;
    double logp_Z = 0;
    double MaxLoglik;
    int Max_m;
    vector <int > Max_Z;
    vector <bool> Max_B;
    double Max_b, Max_n, Max_c;
    
    vector< double >  trace_loglik;  //P(X|Z, TM, r)
    vector< double >  trace_full_loglik;  //P(X,Z,r|TM)
    vector< vector<int> > trace_Z;  //iteration * Species * {0,1,2}, only for init
    vector< vector<bool> > trace_GB;
    vector <double> trace_n_rate;  //iter*element
    vector <double> trace_c_rate;
    vector<double> trace_n_GB;
    //vector <double> trace_c_GB;
    
    int accept_n_rate = 0;
    int accept_c_rate = 0;  //how many accepted in current cycle
    int accept_n_GB;
    int accept_c_GB = 0;
    
    double prop_n;  //for adaptive MCMC, changed by acceptance rate
    double prop_c;
    double prop_nGB;  //for adaptive MCMC, changed by acceptance rate
    vector<double> prop_nGB2;
    double prop_cGB;
    
    
    
    // GSL random number generator
    gsl_rng * RNG;
    
    time_t last_time;
    
    unsigned long int seed;
    
    bool verbose; //= false;

    
public:
    
   
    
    BPP_C(int c, PhyloProf _prof, BPP& bpp, char gapchar, double missing_thres, bool & filter, bool _verbose)//, double _indel)
    {
        
        RNG = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(RNG, bpp.seed);
        
        num_burn = bpp.num_burn;   // num of burn-in updates
        num_mcmc = bpp.num_mcmc;   // num of MCMC updates
        adaptive_freq = bpp.num_thin;   // num of updates between two samples, 100
        
       
        N = bpp.N;
        CC = c;
        GG =bpp.element_size[c];
        //Tg = vector<vector<int> > (GG, vector<int>(N, -1));
        root = N-1;
        S = bpp.S;
        
        lambda = vector< vector<vec> >(GG, vector<vec>(N, zeros<vec>(bpp.num_base)));
        message = vector< vector<vec> >(GG, vector<vec>(N, zeros<vec>(bpp.num_base)));
        //ambiguousS_null = vector<vector<vec> > (GG, vector<vec>(S,zeros<vec>(bpp.num_base)));
        
        
        ratio0 = bpp.ratio0;  //no use
        ratio1 = bpp.ratio1; //no use
        
        verbose = _verbose;
        
        posterior_Tg = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));  //variational posterior prob of I(Tp = a, Ts = b); for root, only use column 0
        posterior_Tg_next = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));  // for q(t+1 -> t+2) at t
        
        log_cache_TM_neut = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_cons = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_null = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        
        log_cache_TM_neut_B = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_cons_B = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_null_B = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        
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
                        lambda[g][s][0] = 0; break; //Tg[g][s] = 0;
                    case 'c':
                        lambda[g][s][1] = 0; break; //Tg[g][s] = 1;
                    case 'g':
                        lambda[g][s][2] = 0; break; //Tg[g][s] = 2;
                    case 't':
                        lambda[g][s][3] = 0; break;//Tg[g][s] = 3; 
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
                        break; //Tg[g][s] = 4; break;
                    default:
                        for(int b =0;b<bpp.num_base;b++) lambda[g][s][b] = 0;
                        //Tg[g][s] = bpp.num_base;
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
        
        
        log_pi = bpp.log_pi;
        
  
        // from initMCMC
        trace_loglik = vector<double >(num_burn+num_mcmc, 0); //P(X|Z, r)
        trace_full_loglik = vector<double >(num_burn+num_mcmc, 0); //P(X, Z, r)
        trace_Z = vector<vector<int> >(num_burn+num_mcmc, vector<int>(N,1));
        trace_GB = vector<vector<bool> >(num_burn+num_mcmc, vector<bool>(N,0));
        trace_n_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_c_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_n_GB = vector<double> (num_burn+num_mcmc, 0.0);
      //  trace_c_GB = vector<double >(num_burn+num_mcmc, 0);

        //getEmission_ambig();
        
        
        prop_c = (double) 1.2; //1.1 0.01; //1/GG;
        prop_n = (double) 1.2; //1.2 0.01; //1/GG ;//* 10;
        //prop_cGB = (double) 1/GG ;
        prop_nGB = 0.2; //1.2 0.4 for real case; //vector<double>(N, 0.1); //(double)1/GG);
        //prop_nGB2 = vector<double>(N, 0.5);  // prob to transit from 1->0
        accept_n_GB = 0; //vector<int>(N, 0);
        
        
        for(int s=0; s<N; s++)
        {
            //log_cache_TM_null[s] = bpp.log_cache_TM_null[s]; // already initialized in hpp, will do in below because change of distance
            log_TM_Int[s] = bpp.log_TM_Int[s];
        }
        
       
        for(vector<int>:: iterator it = upper_c.begin(); it!=upper_c.end();it++)
        {
            log_TM_Int[*it](1,1) = 0;
            log_TM_Int[*it](2,1) = log(0);
        }

    }
    
    
    
    
    ~BPP_C()
    {
        gsl_rng_free(RNG);
    }
    mat getRate(mat Q0, double s, double b);
    mat getTransition(int s, double sel, double propos, BPP & bpp);
    
    void getSubtree(int root, set<int>& child, vector<int> & visited_init);
    void getSubtree(int root, vector<int> & visited_init);
    void getSubtree_missing(int root, set<int>& upper, int child);
    
    //void getEmission_ambig();
    void initMCMC(int iter, BPP&bpp, int resZ);
    void Update_Tg(int g, vector<bool> visited,BPP& bpp, bool tosample);
    void getUpdateNode(vector<int> changedZ, vector<bool> & visited_init);
    void MonitorChain(int m, bool move_r, int& accept,int& accept2, BPP &bpp, const double add_loglik, const int resZ) ;
    void getUpdateNode(bool neut, bool bg, vector<bool> & visited_init);
    double ELBO();
    double log_f_Xz(vector<bool> visited, vector<vector<vec>>& lambda_tmp, vector<vector<vec>>& message_tmp, bool neut,bool gBC,  double propos, BPP& bpp);
    void update_cache_TM(vector<mat> & Qs, vector<vector<mat >>& log_cache_TM_tmp);
    void init_cache_TM(vector<mat> & Qs);
    //double sample_rate(int resZ, double old_rate, bool neut, vector<bool> visited, double & loglik_old, BPP& bpp, int M =1, bool adaptive = true, double adaptive_factor = 1);
    bool sample_rates(BPP& bpp,double & loglik_old, int resZ, int M =1);
    bool update_rates(BPP& bpp, double & elbo_old, int M = 1, bool adaptive = true, double adaptive_factor = 0.2);
    
   // double sample_gBC(double old_rate, vector<bool> visited, double & loglik_old, BPP& bpp, int M =1, bool adaptive = true, double adaptive_factor = 1);
    //double sample_gBC0(int s, double old_rate, BPP& bpp,  int M =1, bool adaptive = true, double adaptive_factor = 1);
    void Gibbs(int iter, BPP &bpp, ofstream & outZ, string output_path,string output_path2,int resZ, bool UpR, bool UpHyper, double lrate_prop, double grate_prop);
  //  vector<int> Update_ZB_subtree(BPP& bpp, bool tosample);
    void integrate_ZB(BPP &bpp, bool gb, bool neut, vector<mat>& Qs, vector<vector<double>> & log_emission_tmp, vector<vector<mat>> & log_cache_TM_tmp);
    vector<int> Update_ZB();

    void Output_init(string output_path, string output_path2, BPP& bpp,ofstream& outZ, int resZ);
    void Output_sampling(int iter, string output_path2, BPP &bpp, int resZ);
    void getEmission_update(bool neut, bool gb, vector<vector<double>> & log_emission_tmp,  vector<vector<mat >>& log_cache_TM_tmp, BPP &bpp, int restore);
    double initial_prob_back(BPP &bpp);
    
    vector<int>  Move_Z(int & propConf, int & revConf,  int & changeZ);
    //double log_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut,  vector<mat> & log_cache_TM_cons);
    double Update_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut, vector<mat> & log_cache_TM_cons, vector<bool> & visited);
    void log_f_Z(vector<int>& Z, vector<mat> & log_Int, double & MH_ratio_g, double & MH_ratio_l);
    
    //void Eval(BPP&bpp,int resZ); //, int numH,int numHZ);
    set<int> getShifts(vector<bool> & Z); 
    void Eval2(BPP&bpp, int resZ);
    double prior_Z(vector<int> & tmpZ) ;
    double prior_Z_subtree(vector<int> & tmpZ) ;
    void adaptive_MCMC(double adaptive_factor, int M =1);
    
};


#endif /* bpp_c_hpp */
