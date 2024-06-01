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
#include <random>

#include "newick.h"
#include "profile.h"
#include "bpp.hpp"
#include "utils.h"
#include "genetree.hpp"


using namespace std;
using namespace arma;

struct simpletree
{
    int count;                       // number of occurence
    vector< double > distances;      // average branch distances
    vector<vector<int>> children_gene;
    vector<int> parent_gene;
    int root;
    vector<string> node_names;     //Han*: name internal node by its descendents in newick format
    void printTree(int s, BPP & bpp, std::stringstream & buffer)
    {
        if (children_gene[s][0] == -1)
        {
            buffer << bpp.species_names[s] << ":"<< distances[s];
        }
        else
        {
            buffer << "(";
            for(int i =0;i <2; i++)
            {
                int child = children_gene[s][i];
                
                printTree(child, bpp, buffer);
                if(i==0) buffer << ",";
                
            }
            
            if(parent_gene[s] < bpp.N)
            {
                buffer << "):" << distances[s];
            }else{
                buffer << ");" ;
            }
        }
    }
};

class BPP_C
{
    private:

    int CC;  //current number of elements
    int GG_block;

    int S;
    int N;

    //vector< vector<vec > > ambiguousS_null;  //base * species# * 4;

    int    (*children2)[2];  //children from pruned tree
    int root;
    //Han*: change to vector
    vec pi;
    vec log_pi;
    //vector<double> pi;
    //vector<double> log_pi;
    vector<double> prior_dir_param;
    //Han*: Q a-f:
    vec inst_rate;
    mat eigenvec; vec eigenval;
    mat eigeninv;
    mat submat;

    
    vector< vector<vec>>  lambda;  //GG * S * vec(4)
    vector < vector<int>> Tg; //current history during each speciation, GG*S
    

    vector<int> nodes;
    //vector<int> internal_nodes;
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
    int adaptive_freq; //= 100;

    vector<int> Z ; //N *0, accelerate(time,0-1), loss(-1)
    vector<int> fixZ ;
   
    //buffer P(zt|zt-1,Hg)
    vector<vec> log_prob_back;

    //get P(X|theta)
    vector< vector<double> >  log_emission; //( N-1)*2

    // MCMC updating states
    int m;                                      // current MCMC step

    vector<mat> log_TM_Int;

    vec prior_z;
    
    // samples to output
    double MaxLoglik;
    int Max_m;
    vector <int > Max_Z;
    string Max_GT;
    //Han*: save running max for pi and acgt counts (for posterior parameters)
    vector<double> Max_pi;

    vector< double >  trace_loglik;  //P(X|Z, TM, r)
    vector< double >  trace_full_loglik;  //P(X,Z,r|TM)
    vector< vector<int> > trace_Z;  //iteration * Species * {0,1,2}, only for init
    vector <double> trace_n_rate;  //iter*element
    vector <double> trace_c_rate;
    vector<double> trace_logNormratio;
    map<string, simpletree> trace_genetree; // record how many time and average branch length

    vector <double> trace_l_rate;  //iter*element
    vector <double> trace_l2_rate;
    vector <double> trace_g_rate;

    //Han*: save MCMC traces
    vector<vector<double>> trace_pi;
    vector<double> trace_indicator;
    vector<int> trace_GTtopChg;

    int accept_n_rate = 0;
    int accept_c_rate = 0;  //how many accepted in current cycle

    double prop_n;  //for adaptive MCMC, changed by acceptance rate
    double prop_c;

    double consToMis;
    double nconsToMis;

    // GSL random number generator
    gsl_rng * RNG;

    time_t last_time;

    unsigned long int seed;
    unsigned long int seed2;
    GTree* gtree;
    friend class GTree;

public:
    int GG; //base pairs current elements
    bool failure = false;
    bool verbose;
    bool verboseGT;
    int idblk_count=0; //length of 1st blck of all identical bp across sp.
    
    BPP_C(int c, PhyloProf _prof, BPP& bpp, char gapchar, double missing_thres, bool & filter, bool _verbose, bool _verboseGT, double _consToMis, int blocks = 20, bool prune=0, double revgap=0, int min_length =50, double _nconsToMis = 0.5)//, double _indel)
    {

        RNG = bpp.RNG;

        num_burn = bpp.num_burn*bpp.num_thin;   // num of burn-in updates
        num_mcmc = bpp.num_mcmc*bpp.num_thin;   // num of MCMC updates
        if(bpp.num_thin == 1)  // 0517: ???
        {
            adaptive_freq = 100;   // num of updates between two samples, 100
        }else{
            adaptive_freq = 50;
        }
        num_thin = bpp.num_thin;

        consToMis = _consToMis;
        nconsToMis = _nconsToMis;

        N = bpp.N;
        CC = c;
        GG =bpp.element_size[c];
        root = N-1;   // for gene tree
        S = bpp.S;
        
        children2    = new int[N][2];
        for(int i=0; i<N; i++)
        {
            children2[i][0] = bpp.children[i][0];
            children2[i][1] = bpp.children[i][1];
        }
        
        gtree = new GTree(N, GG, S, RNG);

        if(GG < min_length)
        {
            filter = true;
            return;
        }
        
        GG_block = blocks; //ceil((double)GG/blocks);
        int tot = ceil((double)(GG - 15)/GG_block);
        // if(tot>5) tot=5;
        // GG_block=ceil((double)GG/(tot));

        //cout << CC << ": " << GG << " bp" << endl;

        lambda = vector< vector<vec> > (GG, vector<vec>(N));
        Tg = vector < vector< int> > (GG, vector<int>(N));
        

        ratio0 = bpp.ratio0;  //no use
        ratio1 = bpp.ratio1; //no use

        verbose = _verbose;
        verboseGT = _verboseGT;

        log_TM_Int = vector<mat >(N, zeros<mat>(3,3));

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

        //Han*:
        prior_dir_param=bpp.prior_dir_param;        

        if(verbose) cout << "Init lambda" <<endl;
        vector<int> num_missing = vector<int> (S,0);

        //get leaves lambda
        int st = bpp.element_start[CC];
        // shuffle sequence
        vector<int> myvector;
        for (int i=0; i<GG; ++i) myvector.push_back(i);
        std::shuffle ( myvector.begin(), myvector.end(), bpp.twister2);//default_random_engine(bpp.seed)

        vector < vector< int> > Tg2 = vector < vector< int> > (GG, vector<int>(N));
        for(int s=0; s<S; s++){
            string y = _prof.X[s].substr(st, st+GG);
            for(int g=0; g<GG; g++){
                int gs=myvector[g];
                Tg2[g][s] = -1;
                switch (y[gs])
                {
                    case 'a':
                        Tg2[g][s] = 0; break;
                    case 'c':
                        Tg2[g][s] = 1; break;
                    case 'g':
                        Tg2[g][s] = 2; break;
                    case 't':
                        Tg2[g][s] = 3; break;
                    case 'r':
                        break;
                    case 'y':
                        break;
                    case 'k':
                        break;
                    case 'm':
                        break;
                    case 's':
                        break;
                    case 'w':
                        break;
                    default:
                        if(y[gs] == gapchar || y[gs]=='n' || y[gs]=='*') {Tg2[g][s] = 4; break;}
			            Tg2[g][s] = 5;
                }
            }
        }
        
        vector<int>::iterator ip;
        vector<int> tmp_tg;
        vector<int> myvector2={}; //put all identical bp blocks in front.
        for(int g=0; g<GG; g++){
            tmp_tg=Tg2[g];
            std::sort(tmp_tg.begin(),tmp_tg.end());
            ip=std::unique(tmp_tg.begin(),tmp_tg.end());
            tmp_tg.resize(std::distance(tmp_tg.begin(),ip));
            if(tmp_tg.size()==1){
                myvector2.insert(myvector2.begin(),g);
                idblk_count+=1;
            }else if(tmp_tg.size()==2){
                if(tmp_tg== vector<int>{0,4} || tmp_tg==vector<int>{0,5} || tmp_tg==vector<int>{1,4} || 
                tmp_tg==vector<int>{1,5} || tmp_tg==vector<int>{2,4} || tmp_tg==vector<int>{2,5} || 
                tmp_tg==vector<int>{3,4} || tmp_tg==vector<int>{3,5} || tmp_tg==vector<int>{4,5}){
                    myvector2.insert(myvector2.begin(),g);
                    idblk_count+=1;
                }else{
                    myvector2.push_back(g);
                }
            }else if(tmp_tg==vector<int> {0,4,5} ||tmp_tg==vector<int> {1,4,5} ||tmp_tg==vector<int> {2,4,5} ||tmp_tg==vector<int> {3,4,5}){
                    myvector2.insert(myvector2.begin(),g);
                    idblk_count+=1;                
            }else{
                myvector2.push_back(g);
            }
        }
        //cout<<"number of idblk is "<<idblk_count<<". myvec2 length "<<myvector2.size()<<". GG="<<GG<<endl;

        for(int s=0; s<S; s++)
        {
            //const char* y = bpp.X[c][s].c_str();
            string y = _prof.X[s].substr(st, st+GG);

            for(int g=0; g<GG; g++) {
                int gs = myvector2[g];
                lambda[g][s] = zeros<vec>(bpp.num_base);
                lambda[g][s].fill(LOG_ZERO);
                Tg[g][s] = -1;
                switch (y[gs])
                {
                    case 'a':
                        lambda[g][s][0] = 0; Tg[g][s] = 0; break;
                    case 'c':
                        lambda[g][s][1] = 0; Tg[g][s] = 1; break;
                    case 'g':
                        lambda[g][s][2] = 0; Tg[g][s] = 2; break;
                    case 't':
                        lambda[g][s][3] = 0; Tg[g][s] = 3; break;
                    case 'r': // a/g
                        lambda[g][s][0] = 0;
                        lambda[g][s][2] = 0;
                        break;
                    case 'y': // c/t
                        lambda[g][s][1] = 0;
                        lambda[g][s][3] = 0;  break;
                    case 'k': // g/t
                        lambda[g][s][2] = 0;
                        lambda[g][s][3] = 0;  break;
                    case 'm': // a/c
                        lambda[g][s][0] = 0;
                        lambda[g][s][1] = 0;  break;
                    case 's': // g/c
                        lambda[g][s][1] = 0;
                        lambda[g][s][2] = 0;  break;
                    case 'w': // a/t
                        lambda[g][s][0] = 0;
                        lambda[g][s][3] = 0;  break;
                    default:
                        for(int b =0;b<bpp.num_base;b++) lambda[g][s][b] = 0;
                        //if(y[gs] == gapchar) {Tg[g][s] = 4; break;}
                        if(y[gs] == gapchar || y[gs]=='n' || y[gs]=='*') {Tg[g][s] = 4; break;}
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
                    if(g<idblk_count) idblk_count = idblk_count -1;
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
        //cout<<"After filter, number of idblk is "<<idblk_count<<". GG="<<GG<<endl;
        if(idblk_count<=5) idblk_count=0; 
        
        gtree->GG = GG; 
        
        // reset GG_block if GG is too long
        //int tot = ceil((double)(GG - 10)/GG_block);
        /*
        if(tot > 5)
        {
            GG_block = ceil((double)(GG - 10)/5); // need blocks >=5!!
        }
        */

        // get gap species after filtering
        for(int s=0; s<S; s++)
        {
            for(int g= 0; g < GG; g++)
            {
                if( Tg[g][s] == 4) //10Nov: >=4?
                {
                    num_missing[s] ++;
                }

            }
        }

        

        // set missing only for extant species and upper
        missing = vector<bool>(N, false);       
        // set missing for all, will be used in getUpdateNode in sample_rate
        for(int s=S; s<N; s++)
        {
            int* p = bpp.children[s];
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
            for(int s=S; s<N; s++)
            {
                int* p = bpp.children[s];
                for(int cc =0 ;cc < 2; cc++)
                {
                    if(missing[p[cc]])
                    {
                        children2[s][cc] = -1;
                    }
                }
            }            
        }else{
            //similar to above, but only cares about nodes on lineages leading to outgroup
            getUppertree_missing(bpp); // get children2: if kid[s] is missing, then s does not have this kid, (reflected in children2)
        }
        
         //init gtree lambda/Tg
        for(int s = 0; s< S; s++)
        {
            gtree->lambda[s][s] = vector<mat>(GG);
            gtree->Tg[s][s] = vector<int>(GG);
            
            for(int g = 0; g < GG; g++)
            {
                gtree->lambda[s][s][g] = zeros<mat>(bpp.num_base, 3);
                gtree->lambda[s][s][g].col(2) = lambda[g][s];
                gtree->Tg[s][s][g] = Tg[g][s];
            }
        }
        
        lambda.clear();
        Tg.clear();
       
        if(verbose) cout << "root: " << root << endl;

        getSubtree(root, nodes);  // will use children2 (only here)
        for(vector<int>::iterator it = nodes.begin(); it < nodes.end(); it++)
        {
            if(bpp.upper.find(*it)!=bpp.upper.end())
            {
                upper_c.push_back(*it);
            }

            if(bpp.upper_conserve.find(*it)!=bpp.upper_conserve.end())
            {
                upper_conserve_c.push_back(*it);
            }

        }

        //Han*
        pi=bpp.pi;
        log_pi = bpp.log_pi;
        inst_rate=bpp.inst_rate;
        submat=bpp.submat;
        eigenvec=bpp.eigenvec;
        eigenval=bpp.eigenval;
        eigeninv=bpp.eigeninv;


        // from initMCMC
        trace_loglik = vector<double >(num_burn+num_mcmc, 0); //P(X|Z, r)
        trace_full_loglik = vector<double >(num_burn+num_mcmc, 0); //P(X, Z, r)
        trace_Z = vector<vector<int> >(num_burn+num_mcmc, vector<int>(N,1));
        trace_n_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_c_rate = vector<double >(num_burn+num_mcmc, 0);
        trace_l_rate = vector<double >(num_burn+num_mcmc, bpp.cur_lrate[CC]);
        trace_g_rate = vector<double >(num_burn+num_mcmc, bpp.cur_grate[CC]);
        trace_l2_rate = vector<double >(num_burn+num_mcmc, bpp.cur_lrate2[CC]);
        trace_logNormratio = vector<double>(num_burn+num_mcmc, 0);


        //Han*:
        trace_pi=vector<vector<double>> (num_burn+num_mcmc,vector<double>(4,0.25)); //all inialized to pi from input
        //Han**: 0322-debug
        trace_indicator=vector<double>(num_burn+num_mcmc,0);
        trace_GTtopChg = vector<int>(num_burn+num_mcmc,0);
        
        //initalize Z, all Z are 1 except the root
        for(int res =0; res < 3; res ++)
        {
            bpp.cur_Z[res][CC] = vector<int>(N,-1);
            for(vector<int>::iterator it = nodes.begin();it<nodes.end();it++) {
                bpp.cur_Z[res][CC][*it] = 1;
            }
            bpp.cur_Z[res][CC][root] = 0; 
        }
        log_emission = vector<vector<double> >(N, vector<double>(3,0));
    }

    ~BPP_C()
    {
        //gsl_rng_free(RNG);
        delete [] children2;
        delete gtree;
    }

    void getSubtree(int root, vector<int> & visited_init);
    void getUppertree_missing(BPP & bpp);
    void initMCMC(int iter, int max_iter, BPP&bpp, int resZ, bool prune, bool fixtree = false); 
    void getUpdateNode(vector<int> changedZ, vector<bool> & visited_init, BPP & bpp);
    void getUpdateNode(bool neut,vector<bool> & visited_init, BPP & bpp);
    void MonitorChain(int m,int iter, int max_iter, BPP &bpp, double ind_prop, const double add_loglik, const int resZ, bool recordtree = true) ;
    double sample_rate(int indictor,int iter, int max_iter, vector<int> lens, int resZ, double old_rate, bool neut, vector<bool> visited, vec & loglik_old, BPP& bpp, int M =1, bool adaptive = true, double adaptive_factor = 0.3);
    void Gibbs(int iter, int max_iter, BPP &bpp, ofstream & outZ, string output_path,string output_path2,int resZ, bool UpR, bool UpHyper, double lrate_prop, double grate_prop, bool WL = true);
    vector<int> Update_Z(int len, BPP &bpp);
    void Output_init(string output_path, string output_path2, BPP& bpp, ofstream& outZ, ofstream& out_tree, int mod_GT);
    void Output_tree(int iter, string outpathG, BPP &bpp, int resZ);
    void Output_sampling(int iter, string output_path2, BPP &bpp, int resZ);
    void getEmission(int len, BPP & bpp);
    vector<int>  Move_Z(int & propConf, int & revConf,  int & changeZ);
    void sample_transition( double  & gr, double  & lr, double  & lr2, BPP& bpp);
    void simulate(BPP& bpp, PhyloProf & profile, char gapchar, bool prune);
    double priorP_Z(BPP & bpp);
    void getQmat(vec & statR, vec & instR);
    int sample_pi(int m, int indictor, vector<int> lens, int resZ, double piA_old, vec &loglik_old, double pi_delta , BPP &bpp, vector<double> obsCount, int resample_burnin);
    void Output_GTsampling(string output_path2,BPP& bpp, int resZ);
    void printSptree(BPP& bpp, int l=0);
};

#endif /* bpp_c_hpp */
