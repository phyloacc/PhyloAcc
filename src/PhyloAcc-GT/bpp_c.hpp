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

#include "../PhyloAcc-common/newick.h"
#include "../PhyloAcc-common/profile.h"
#include "../PhyloAcc-common/bpp_constructor.h"
#include "../PhyloAcc-common/rng.h"
#include "bpp.hpp"
#include "../PhyloAcc-common/utils.h"
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

    int CC = 0;  //current number of elements
    int GG_block = 0;

    int S = 0;
    int N = 0;

    //vector< vector<vec > > ambiguousS_null;  //base * species# * 4;

    int    (*children2)[2] = nullptr;  //children from pruned tree
    int root = -1;
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
    double prior_l_a = 0.0, prior_l_b = 0.0;
    double prior_l2_a = 0.0, prior_l2_b = 0.0;
    double prior_g_a = 0.0, prior_g_b = 0.0;
    double ratio0 = 0.0;
    double ratio1 = 0.0;
    int num_burn = 0;   // num of burn-in updates
    int num_mcmc = 0;   // num of MCMC updates
    int num_thin = 0;   // num of updates between two samples
    int adaptive_freq = 0; //= 100;

    vector<int> Z ; //N *0, accelerate(time,0-1), loss(-1)
    vector<int> fixZ ;
   
    //buffer P(zt|zt-1,Hg)
    vector<vec> log_prob_back;

    //get P(X|theta)
    vector< vector<double> >  log_emission; //( N-1)*2

    // MCMC updating states
    int m = 0;                                      // current MCMC step

    vector<mat> log_TM_Int;

    vec prior_z;
    
    // samples to output
    double MaxLoglik = -INFINITY;
    int Max_m = 0;
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

    double prop_n = 0.0;  //for adaptive MCMC, changed by acceptance rate
    double prop_c = 0.0;

    double consToMis = 0.0;
    double nconsToMis = 0.0;

    // GSL random number generator
    gsl_rng * RNG = nullptr;

    time_t last_time = 0;

    unsigned long int seed = 0;
    unsigned long int seed2 = 0;
    GTree* gtree = nullptr;
    friend class GTree;

public:
    int GG = 0; //base pairs current elements
    bool failure = false;
    bool verbose = false;
    bool verboseGT = false;
    int idblk_count=0; //length of 1st blck of all identical bp across sp.
    
    BPP_C(int c, PhyloProf _prof, BPP& bpp, char gapchar, double missing_thres, bool & filter, bool _verbose, bool _verboseGT, double _consToMis, int blocks = 20, bool prune=0, double revgap=0, int min_length =50, double _nconsToMis = 0.5, int chain_index = 0)//, double _indel)
    {

        RNG = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(RNG, phyloacc::DeriveSeed(bpp.seed, phyloacc::ProgramKind::GT,
                                              chain_index, c, blocks,
                                              phyloacc::RngStream::WorkerGsl));

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
        
        unsigned long gene_tree_seed = phyloacc::DeriveSeed(
            bpp.seed2, phyloacc::ProgramKind::GT, chain_index, c, blocks,
            phyloacc::RngStream::GtGeneTreeShuffle);
        gtree = new GTree(N, GG, S, RNG, phyloacc::MakeTwister(gene_tree_seed));

        if(GG < min_length)
        {
            filter = true;
            return;
        }
        
        GG_block = blocks; //ceil((double)GG/blocks);
        // if(tot>5) tot=5;
        // GG_block=ceil((double)GG/(tot));

        //cout << CC << ": " << GG << " bp" << endl;

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

        //get leaves lambda
        int st = bpp.element_start[CC];
        // shuffle sequence
        vector<int> myvector;
        for (int i=0; i<GG; ++i) myvector.push_back(i);
        unsigned long site_shuffle_seed = phyloacc::DeriveSeed(
            bpp.seed, phyloacc::ProgramKind::GT, chain_index, c, blocks,
            phyloacc::RngStream::GtSiteShuffle);
        std::mt19937 site_shuffle_rng = phyloacc::MakeTwister(site_shuffle_seed);
        std::shuffle ( myvector.begin(), myvector.end(), site_shuffle_rng);//default_random_engine(bpp.seed)

        vector < vector< int> > Tg2 = vector < vector< int> > (GG, vector<int>(N));
        for(int s=0; s<S; s++){
            for(int g=0; g<GG; g++){
                int gs=myvector[g];
                Tg2[g][s] = phyloacc::EncodeLeafState(
                    _prof.X[s][st + gs], gapchar, phyloacc::MissingBasePolicy::GapNStar);
            }
        }
        
        vector<int> tmp_tg;
        vector<int> myvector2={}; //put all identical bp blocks in front.
        for(int g=0; g<GG; g++){
            tmp_tg=Tg2[g];
            if(phyloacc::IsSimpleOrMissingLeafPattern(tmp_tg)){
                    myvector2.insert(myvector2.begin(),g);
                    idblk_count+=1;
            }else{
                myvector2.push_back(g);
            }
        }
        //cout<<"number of idblk is "<<idblk_count<<". myvec2 length "<<myvector2.size()<<". GG="<<GG<<endl;

        phyloacc::LeafEncoding leaf_encoding = phyloacc::EncodeLeafAlignment(
            _prof.X, st, GG, S, N, bpp.num_base, gapchar,
            phyloacc::MissingBasePolicy::GapNStar, myvector2);
        lambda = leaf_encoding.lambda;
        Tg = leaf_encoding.tg;

        phyloacc::ColumnFilterResult column_filter = phyloacc::RemoveHighMissingColumns(
            lambda, Tg, S, revgap, min_length, &idblk_count);
        if(column_filter.filtered)
        {
            filter = true;
            return;
        }

        GG = column_filter.length;
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

        vector<int> num_missing = phyloacc::CountMissingBySpecies(Tg, S);

        

        // set missing only for extant species and upper
        // set missing for all, will be used in getUpdateNode in sample_rate
        missing = phyloacc::BuildMissingNodes(S, N, bpp.children, num_missing, missing_thres, GG);

        if(phyloacc::ConservedMissingExceeds(missing, bpp.conservedgroup, bpp.conserve_prop))
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
        phyloacc::CollectUpperNodesInSubtree(nodes, bpp.upper, bpp.upper_conserve,
                                             upper_c, upper_conserve_c);

        //Han*
        pi=bpp.pi;
        log_pi = bpp.log_pi;
        inst_rate=bpp.inst_rate;
        submat=bpp.submat;
        eigenvec=bpp.eigenvec;
        eigenval=bpp.eigenval;
        eigeninv=bpp.eigeninv;


        // from initMCMC
        phyloacc::BppCTraceBuffers traces = phyloacc::InitializeBppCTraceBuffers(
            num_burn+num_mcmc, N, bpp.cur_lrate[CC], bpp.cur_lrate2[CC], bpp.cur_grate[CC]);
        trace_loglik = traces.trace_loglik; //P(X|Z, r)
        trace_full_loglik = traces.trace_full_loglik; //P(X, Z, r)
        trace_Z = traces.trace_z;
        trace_n_rate = traces.trace_n_rate;
        trace_c_rate = traces.trace_c_rate;
        trace_l_rate = traces.trace_l_rate;
        trace_g_rate = traces.trace_g_rate;
        trace_l2_rate = traces.trace_l2_rate;
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
        log_emission = traces.log_emission;
    }

    ~BPP_C()
    {
        delete [] children2;
        delete gtree;
        gsl_rng_free(RNG);
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
