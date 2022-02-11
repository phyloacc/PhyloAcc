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
    vector<int> upper_conserve_c;
    
   // double prior_glr[3];
    double prior_l_a, prior_l_b;
     double prior_l2_a, prior_l2_b;
     double prior_g_a, prior_g_b;
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
    
    vec prior_z;
    
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
    
    vector <double> trace_l_rate;  //iter*element
    vector <double> trace_l2_rate;
    vector <double> trace_g_rate;
    
    int accept_n_rate = 0;
    int accept_c_rate = 0;  //how many accepted in current cycle
    
    double prop_n;  //for adaptive MCMC, changed by acceptance rate
    double prop_c;
    
    
    
    // GSL random number generator
    gsl_rng * RNG;
    
    time_t last_time;
    
    unsigned long int seed;

    
public:
    bool failure = false;
    bool verbose; 
    
    BPP_C(int c, PhyloProf _prof, BPP& bpp, char gapchar, double missing_thres, bool & filter, bool _verbose, double consToMis, bool prune=0, double revgap=0, int min_length =50, double nconsToMis = 1)//, double _indel)
    {
        
        RNG = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(RNG, bpp.seed);
        
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
        
        Tg = vector<vector<int> > (GG, vector<int>(N, -1));
        root = N-1;
        S = bpp.S;
        
        lambda = vector< vector<vec> >(GG, vector<vec>(N, zeros<vec>(bpp.num_base)));
        
        
        ratio0 = bpp.ratio0;  //no use
        ratio1 = bpp.ratio1; //no use
        
        verbose = _verbose;

        
        log_cache_TM_neut = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_cons = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        log_cache_TM_null = vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base));
        
        log_TM_Int = vector<mat >(N, zeros<mat>(3,2));
        
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
        vector<int> num_missing = vector<int> (S,0);
        
        int st = bpp.element_start[CC];
        //get leaves lambda
        for(int s=0; s<S; s++)
        {
            //const char* y = bpp.X[c][s].c_str();
            string y = _prof.X[s].substr(st, st+GG);
            
            for(int g=0; g<GG; g++){
                lambda[g][s].fill(LOG_ZERO);
//                if( y[g]== gapchar || y[g]=='n')
//                {
//                    num_missing[s] ++;
//                    //Tg[g][s] = 4;  // for missing base pair
//                }
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
                    // case gapchar:
                    //     if(bpp.num_base <= 4){
                    //         for(int b =0;b<bpp.num_base;b++) lambda[g][s][b] = 0;
                    //     }
                    //     else{
                    //         lambda[g][s][4] = 0;
                    //     }
                    //     Tg[g][s] = 4; break;
                    default:
                        for(int b =0;b<bpp.num_base;b++) lambda[g][s][b] = 0;
                        if(y[g] == gapchar) {Tg[g][s] = 4; break;}
                        Tg[g][s] = 5; 
                }
            }
            
        }
        
        
        // remove columns with nearly all gaps
        if(revgap < 1)
        {
            //remove bases with 'n'/'*' or gap in more than 80% species
            vector<int> missingBase;
            for(int g=0; g<GG; g++)
            {
                int mis = 0;
                for(int s=0; s<S; s++){
                    if(Tg[g][s] >= 4)
                    {
                        mis++;
                    }
                }
                if(mis > S*revgap)
                {
                    missingBase.push_back(g);
                }
                
            }
            
            if(GG - missingBase.size() < min_length)
            {
                filter = true;
                return;
            }
            
            for(vector<int>::reverse_iterator it = missingBase.rbegin(); it!= missingBase.rend(); it ++)
            {
                lambda.erase(lambda.begin() + *it);
                Tg.erase(Tg.begin() + *it);
            }
            
            GG = lambda.size();
        }
        
        // get gap species after filtering
        for(int s=0; s<S; s++)
        {
            for(int g= 0; g < GG; g++)
            {
                if( Tg[g][s] == 4)
                {
                    num_missing[s] ++;
                }
                
            }
        }

        ambiguousS_null = vector<vector<vec> > (GG, vector<vec>(S,zeros<vec>(bpp.num_base)));
        
        
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
            
            if(bpp.upper.find(*it)!=bpp.upper.end())
            {
                    upper_c.push_back(*it);
            }
            
            if(bpp.upper_conserve.find(*it)!=bpp.upper_conserve.end())
            {
                upper_conserve_c.push_back(*it);
            }

        }
        
        
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
        trace_loglik = vector<double >(num_burn+num_mcmc, 0); //P(X|Z, r)
        trace_full_loglik = vector<double >(num_burn+num_mcmc, 0); //P(X, Z, r)
        trace_Z = vector<vector<int> >(num_burn+num_mcmc, vector<int>(N,1));
        trace_n_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_c_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_l_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_g_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_l2_rate = vector<double >(num_burn+num_mcmc, 0);

        //getEmission_ambig();
        
        
        
        
        log_emission = vector<vector<double> >(N, vector<double>(3,0));
        
        
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
