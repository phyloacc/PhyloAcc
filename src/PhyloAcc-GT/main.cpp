//
//  main.cpp
//  PhyloAcc
//
//  Created by hzr on 3/8/16.
//  Copyright © 2016 hzr. All rights reserved.
//
#include <dirent.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iomanip>
//#include <omp.h> //libiomp/
#include <string>
#include <armadillo>
#include "profile.h"
#include "newick.h"
#include "newick2.h"
#include "bpp.hpp"
#include "bpp_c.hpp"
#include "utils.h"
#include <math.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace arma;

// parameters, input and output files paths
string params_path;
string phytree_path;
string align_path;
string output_path="";
string output_path2="";
string segment_path;
string id_path="";
string result_prefix="test";
string tree_coal_unit="";

//string refspecies;
string outgroup;
string targetspecies;
string conservegroup; // can't be missing in more than 50%
string deepcoal_species = ""; //species in which deep coal can happen
double conserve_prop = 0.8;

int num_thread = 1;

// running parameters
int num_burn = 200;         // num of burn-in updates is num_burn * num_thin
int num_mcmc = 800;         // num of MCMC updates, num_burn * num_mcmc
int num_thin = 1;    // num of updates between two samples, adaptive frequency = 50 * num_thin
int num_chain; // outer loop of updates Q matrix and hyperparameter of substitution rates


double prep_lrate = 0.5;
double prep_lrate2 = 0.1; //
double prep_grate = 0.8; // initalization

double prior_lrate2_a = 1,prior_lrate2_b = 1 ; // beta prior for lrate2, 0.5
double prior_lrate_a = 1 ,prior_lrate_b = 1 ; // beta prior for lrate, 1,9
double prior_grate_a = 1,prior_grate_b = 1; // beta prior for grate, 3,1


double ratio0 = 0.5; //initial conserved rate,0.5
double ratio1 = 1; // initial accelerated rate
double missing_thres = 0.8;

double nprior_a = 10, nprior_b = 0.2;  //around 1
double cprior_a = 5, cprior_b = 0.04;  //around ratio
int ropt = 1;
double cub = 1;
double nlb = 0.6;

int batch = -1 ;
int seed = 1;
int seed2 = 1;
double indel = 0;
double indel2;
bool sample_indel = 0;
bool sample_hyper = false;
char gapchar = '-';
bool verbose = 0;
double consToMis = 0.5; //for simulation, 0.01;
int block = 15; // 25; //90;
bool prune=0;
double revgap=0.9;  // 1
int min_length = 50;
bool WL = true;
bool simulate=false;
bool verboseGT = 1; //Han*: output trace_genetrees
double br_sample_cutoff = 10.0;
double theta_cutoff = 1.0;

//Han*: add Dirichlet prior for stationary distribution pi:
//vector<double> prior_dir_par(4,10);
vector<double> prior_dir_par(2,10);  //Han: it's beta prior for pi. Keep the old name.

// load the program parameters
void LoadParams(int argc, char* argv[])
{
    cout << "Loading input data and running parameters......" << endl;

    if (argc > 1)
        params_path = string(argv[1]);
    else
        params_path = "/Users/zhu/PhyloAcc/Yan-etal-2022/param_CaseB.txt"; //params2
        //cerr << " No parameter file specified. Exiting." << endl;
        //exit(1);

    cout << "Loading program configurations from " << params_path << "......" <<endl;

    const int BUFF_SIZE = 1024;
    char line_buff[BUFF_SIZE];

    ifstream in_params(params_path.c_str());
    if (!in_params)
    {
        cerr << "Cannot open the parameters file: " << params_path.c_str() << endl;
        exit(1);
    }
    while(!in_params.eof())
    {
        in_params.getline(line_buff, BUFF_SIZE);
        istringstream line_stream(line_buff);
        string tmp; line_stream >> tmp;
        
        // simulate or inference
        if (tmp=="SIMULATE")
            line_stream >> simulate;
        // input and output file paths
        else if (tmp=="PHYTREE_FILE")
            line_stream >> phytree_path;
        else if (tmp=="TREE_IN_COALESCENT_UNIT")
            line_stream >> tree_coal_unit;
        else if (tmp=="ALIGN_FILE")
            line_stream >> align_path;
        else if (tmp=="SEG_FILE")
            line_stream >> segment_path;
        else if (tmp=="ID_FILE")
            line_stream >> id_path;
        else if (tmp=="BATCH")
            line_stream >> batch;
        else if (tmp=="RESULT_FOLDER")
            line_stream >> output_path;
        else if (tmp=="PREFIX")
            line_stream >> result_prefix;
        else if (tmp=="SEED")
            line_stream >> seed;
        else if (tmp=="SEEDS")
            line_stream >> seed2;
        else if (tmp=="INIT_CONSERVE_RATE")
            line_stream >> ratio0;
        else if (tmp=="INIT_ACCE_RATE")
            line_stream >> ratio1;
        else if (tmp=="CONSERVE_PRIOR_A")
            line_stream >> cprior_a;
        else if (tmp=="CONSERVE_PRIOR_B")
            line_stream >> cprior_b;
        else if (tmp=="ACCE_PRIOR_A")
            line_stream >> nprior_a;
        else if (tmp=="ACCE_PRIOR_B")
            line_stream >> nprior_b;
        else if (tmp=="ROPT")
            line_stream >> ropt;
        else if (tmp=="CUB")
            line_stream >> cub;
        else if (tmp=="NLB")
            line_stream >> nlb;

        // running parameters
        else if (tmp=="BURNIN")
            line_stream >> num_burn;
        else if (tmp=="MCMC")
            line_stream >> num_mcmc;
        else if (tmp=="THIN") //ADAPT_FREQ
            line_stream >> num_thin;
        else if (tmp=="INIT_LRATE")
            line_stream >> prep_lrate;
	else if (tmp=="INIT_LRATE2")
            line_stream >> prep_lrate2;
        else if (tmp=="INIT_GRATE")
            line_stream >> prep_grate;
        else if (tmp=="HYPER_LRATE_A")
            line_stream >> prior_lrate_a;
        else if (tmp=="HYPER_LRATE_B")
            line_stream >> prior_lrate_b;
        else if (tmp=="HYPER_LRATE2_A")
            line_stream >> prior_lrate2_a;
        else if (tmp=="HYPER_LRATE2_B")
            line_stream >> prior_lrate2_b;
        else if (tmp=="HYPER_GRATE_A")
            line_stream >> prior_grate_a;
        else if (tmp=="HYPER_GRATE_B")
            line_stream >> prior_grate_b;
        else if (tmp=="CHAIN")
            line_stream >> num_chain;
        else if (tmp=="WL")
            line_stream >> WL;
        else if (tmp=="BLK_WL")
            line_stream >> block;   

        // constraint
        else if (tmp == "OUTGROUP")
            line_stream >> outgroup;
        else if (tmp == "TARGETSPECIES")
            line_stream >> targetspecies;
        else if (tmp == "CONSERVE")
            line_stream >> conservegroup;
        else if (tmp == "CONSERVE_PROP")
            line_stream >> conserve_prop;
        else if (tmp == "GAP_PROP")
            line_stream >> missing_thres;
        //else if (tmp == "REF")
         //   line_stream >> refspecies;
        else if (tmp == "CONSTOMIS")
            line_stream >> consToMis;
        else if (tmp=="BR_SAMPLE_THRESHOLD")
            line_stream >> br_sample_cutoff;
        else if (tmp=="THETA_CUTOFF")
            line_stream >> theta_cutoff;
        else if( tmp=="DEEP_COAL_BRANCH")
            line_stream >> deepcoal_species;

        // treat indel as additional character
        else if (tmp == "GAPCHAR")
            line_stream >> gapchar;
        else if (tmp == "PRUNE_TREE")
            line_stream >> prune;
        else if (tmp == "TRIM_GAP_PERCENT")
            line_stream >> revgap;
        else if (tmp == "MIN_LEN")
            line_stream >> min_length;
        else if (tmp == "INDEL") // not used
            line_stream >> indel;
        else if (tmp == "INDEL2") // not used
            line_stream >> indel2;
        else if(tmp == "SAMPLE_INDEL")  // not used
            line_stream >> sample_indel;
        else if(tmp == "SAMPLE_HYPER")
            line_stream >> sample_hyper;
        else if(tmp == "VERBOSE")
            line_stream >> verbose;
        else if(tmp == "VERBOSE_GENETREE") //Han*: output sampled gene trees
            line_stream>>verboseGT;
        else if(tmp == "NUM_THREAD")
            line_stream >> num_thread;
        else if(tmp != "")
            cout << "Unknown parameter: " << tmp <<endl;

    }
    
    if(prior_lrate2_a==0) prep_lrate2 = 0;

    // trimming file names
    phytree_path = strutils::trim(phytree_path, " \"\t\n");
    align_path = strutils::trim(align_path, " \"\t\n");
    output_path  = strutils::trim(output_path,  " \"\t\n");
    segment_path = strutils::trim(segment_path, " \"\t\n");
    tree_coal_unit = strutils::trim(tree_coal_unit, " \"\t\n");
}

bool DirectoryExists( string pzPath )
{
    if ( pzPath == "") return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath.c_str());

    if (pDir != NULL)
    {
        bExists = true;
        (void) closedir (pDir);
    }

    return bExists;
}

void DispParams(PhyloProf profile, int seed)
{
    double mean_seg_size = 0;
    for(unsigned int c=0; c<profile.C; c++)
        mean_seg_size += (double)(profile.element_pos[c][1] - profile.element_pos[c][0]) / profile.C;
    cout << "  # total length = " << profile.G << " (" << profile.C << ")" << ". # Species = " << profile.S << ". # elements = " << profile.C << ". Mean gene set size = " << mean_seg_size << "." << endl;
    cout << "# Burn-ins = " << num_burn*num_thin << ". # MCMC Updates = " << num_mcmc*num_thin << ". # thin = " << num_thin << ".  RND SEED = " << seed << "." << endl ; //
    cout << "# Threads = " << num_thread << endl << endl;
}

int main(int argc, char* argv[])
{
    time_t start = time(NULL);

    cout << std::fixed << setprecision(4);
    srand(time(NULL));

    // load the program parameters
    LoadParams(argc, argv);

    // check output path
    if(! DirectoryExists(output_path))
    {
    	cout << "output path doesn't exist or empty!" << endl;
    	return 1;
    }
    
    // load the phylogenetic profile
    PhyloProf profile = LoadPhyloProfiles(align_path,segment_path);
    // init and display the running parameters
    DispParams(profile, seed);

    // load the phylogenetic tree
    PhyloTree phytree = LoadPhyloTree(phytree_path); //Han: .subs_rate contains Q
    
    PhyloTree_theta tree2;
    if(tree_coal_unit !=""){
        tree2 = LoadPhyloTree_theta(tree_coal_unit);
        //get thetas
        double theta_cum=0;
        int count_cum=0;
        int N = phytree.nodes_names.size();
        vector<int> pos_cum=vector<int> (N,0);
        for(int i=0; i<(N-1); i++){
            if(phytree.nodes_names[i]!=tree2.nodes_names[i]){
                cout<<"i="<<i<<". tree1_name="<<phytree.nodes_names[i]<<", tree2_name="<<tree2.nodes_names[i]<<endl;
                cerr<<"two trees do not have the same topology"<<endl;
                exit(1);
            }else{
                if(i<phytree.S){
                    phytree.thetas[i]=0;
                }else{
                    //cout<<"\ni="<<i<<": ";
                    if((tree2.distances[i]!=1.0) && (tree2.distances[i]!= 7.0) && (tree2.distances[i]!= 0)){
                        phytree.thetas[i]=2*phytree.distances[i]/tree2.distances[i];
                        if(phytree.thetas[i]>= theta_cutoff){
                            pos_cum[i]=1;
                        }else{
                            theta_cum+=phytree.thetas[i];
                            count_cum+=1;
                        }
                    }else{
                        pos_cum[i]=1;
                    }
                }
            }
        }
        double theta_aver=theta_cum/count_cum;
        pos_cum[N-1]=1;
        for(int i=phytree.S; i<N; i++){
            if(pos_cum[i]==1) phytree.thetas[i]=theta_aver;
            //cout<<"node "<<phytree.nodes_names[i]<<" theta="<<phytree.thetas[i]<<endl;
        }
    }else{
        cout<<"please input a phylogengy with branch length in coalescent unit."<<endl;
        return 1;
    } 


    // create and init the BPP object
    //BPP bpp(0, profile, phytree, output_path, targetspecies, outgroup, conserve_prop, conservegroup, ratio0, ratio1, ropt, cub, nlb, nprior_a, nprior_b, cprior_a, cprior_b, seed, seed2, prep_grate, prep_lrate, prep_lrate2, prior_grate_a, prior_grate_b,prior_lrate_a, prior_lrate_b,prior_lrate2_a, prior_lrate2_b,  indel, indel2, missing_thres, sample_indel,prior_dir_par, br_sample_cutoff);
    BPP bpp(0, profile, phytree, output_path, targetspecies, outgroup, conserve_prop, conservegroup, ratio0, ratio1, ropt, cub, nlb, nprior_a, nprior_b, cprior_a, cprior_b, seed, seed2, prep_grate, prep_lrate, prep_lrate2, prior_grate_a, prior_grate_b,prior_lrate_a, prior_lrate_b,prior_lrate2_a, prior_lrate2_b,  indel, indel2, missing_thres, sample_indel,prior_dir_par, br_sample_cutoff, deepcoal_species);

    //initialize the MCMC sampling
    bpp.InitMCMC(num_burn, num_mcmc, num_thin);
    
    output_path = output_path + "/" + result_prefix ;
    output_path2 = output_path;
    string outpath_Z0 = output_path + "_rate_postZ_M" +to_string(0) +".txt";
    string outpath_Z1 = output_path + "_rate_postZ_M" +to_string(1) +".txt";
    string outpath_Z2 = output_path + "_rate_postZ_M" +to_string(2) +".txt";

    string outpath_hyper = output_path+"_hyper.txt";
    ofstream out_hyper(outpath_hyper.c_str());
    out_hyper << "iter\tnprior_a\tnprior_b\tcprior_a\tcprior_b\tprior_l_a\tprior_l_b\tprior_g_a\tprior_g_b\n";
    out_hyper << 0 << "\t"<< nprior_a<< "\t"<< nprior_b <<"\t"<< cprior_a << "\t"<< cprior_b << "\t"<< prior_lrate_a << "\t"<< prior_lrate_b << "\t"<< prior_grate_a << "\t"<< prior_grate_b <<endl;

    ofstream out_lik;
    if(sample_hyper) {
        string outpath_elem = output_path+ "_elem_lik.txt";
        out_lik.open(outpath_elem.c_str());
        out_lik.precision(8);
        out_lik << "No.\tID\tloglik_Full\tloglik_Max"<<endl;
    }

    ofstream out_Z0(outpath_Z0.c_str());
    ofstream out_Z1(outpath_Z1.c_str());
    ofstream out_Z2(outpath_Z2.c_str());

    // output species name
    string species_name = output_path+"_species_names.txt";
    ofstream out_species(species_name.c_str());

    out_Z0 << "No.\tn_rate\tc_rate\tg_rate\tl_rate\tl2_rate"; out_Z1 << "No.\tn_rate\tc_rate\tg_rate\tl_rate\tl2_rate"; out_Z2 << "No.\tn_rate\tc_rate\tg_rate\tl_rate\tl2_rate";
    for(int s=0; s<bpp.N;s++){
         for(int k=0;k<4;k++){
            out_Z0 <<"\t"<<bpp.nodes_names[s]<<"_"<<k;
            out_Z1 <<"\t"<<bpp.nodes_names[s]<<"_"<<k;
            out_Z2 <<"\t"<<bpp.nodes_names[s]<<"_"<<k;
         }
         out_species << bpp.nodes_names[s] << endl;
    }
	out_Z0 <<endl; out_Z1 <<endl; out_Z2 <<endl;

    out_species.close();
    
    // output tree
    outpath_Z0 = output_path + "_tree_M" +to_string(0) +".txt";
    outpath_Z1 = output_path + "_tree_M" +to_string(1) +".txt";
    outpath_Z2 = output_path + "_tree_M" +to_string(2) +".txt";
    
    ofstream out_tree0(outpath_Z0.c_str());
    ofstream out_tree1(outpath_Z1.c_str());
    ofstream out_tree2(outpath_Z2.c_str());
    
    out_tree0 << "No.\tprop\tgenetree\n";
    out_tree1 << "No.\tprop\tgenetree\n";
    out_tree2 << "No.\tprop\tgenetree\n";

    double lrate_prop = 0.5, grate_prop = 0.5;

    vector<int> ids;
    if(id_path=="")
    {
        if(batch==-1)
        {
          for(int c =0;c<500;c++) //bpp.C,
          {
            ids.push_back(c);
          }
        }else{
          int temp = ceil(bpp.C/3);
          for(int c =batch*temp ;c< (batch+1)*temp;c++)
          {
            if(c >= bpp.C) break;
            ids.push_back(c);
          }

        }

    }else{
        ifstream in_params(id_path.c_str());
        if (!in_params)
        {
            cerr << "Cannot open the id file: " << id_path.c_str() << endl;
            exit(1);
        }
        while(!in_params.eof())
        {
            const int BUFF_SIZE = 1024;
            char line_buff[BUFF_SIZE];
            in_params.getline(line_buff, BUFF_SIZE);
            istringstream line_stream(line_buff);
            string tmp; line_stream >> tmp;
            tmp = strutils::trim(tmp);
            if(tmp=="") continue;
            ids.push_back(atoi(tmp.c_str()));

        }

    }

    cout << ids.size() << " elements to be computed" << endl;

    if(sample_hyper)
    {
        for(int iter  =0; iter<num_chain; iter++)
        {
            cout << "Running MCMC chain " << iter +1 << " ..." << endl;
            // Gibbs sampling
            #pragma omp parallel for schedule (guided) num_threads(num_thread)
            for(std::size_t i = 0; i < ids.size(); i++ )
            {
                int c = ids[i];
                bool filter = false;

                try{
                    BPP_C bppc(c, profile, bpp, gapchar, missing_thres, filter, verbose, consToMis, block, prune, revgap, min_length);  // for individual element
                    if(filter) {
                        if(verbose) cerr << "filter: "<< c <<endl;
                        continue;
                    }

                    bppc.initMCMC(0,5,bpp,1,prune);
                    bppc.Gibbs(0, 4, bpp,out_Z2,output_path,output_path2,1, true, sample_hyper, lrate_prop, grate_prop, false);  // Gibbs run to get Z for each element

                    if(bppc.verbose || bppc.failure) bppc.Output_sampling(iter, output_path2, bpp, 2);
                    bppc.Output_init(output_path,output_path2,bpp,out_Z2, out_tree2, bppc.verbose); //sort rates!!

                }catch (exception& e){
                    cout << c << " Standard exception: " << e.what() << endl;
                }
            }
            bpp.sample_hyperparam(iter, ids, out_hyper);
            bpp.Output_init0(profile,out_lik, ids);

        }
    }else if(simulate)
    {
        for(std::size_t i = 0; i < ids.size(); i++ )
        {
            int c = ids[i];
            bool filter = false;
            try{
                // accelerate in target species
                BPP_C bppc(c, profile, bpp, gapchar, missing_thres, filter, verbose, consToMis, block, prune, revgap, min_length);  // for individual element
                bppc.simulate(bpp, profile, gapchar,prune);
            }catch (exception& e){
                cout << c << " Standard exception: " << e.what() << endl;
            }
        }

        //write out simulate sequence
        bpp.Output_simu(profile, output_path, ids.size());        
    }else{
        // Gibbs sampling
        #pragma omp parallel for schedule (guided) num_threads(num_thread)
        for(std::size_t i = 0; i < ids.size(); i++ ) //
        {
            int c = ids[i];
            bool filter = false;
            
            try{
               
                BPP_C bppc(c, profile, bpp, gapchar, missing_thres, filter, verbose, verboseGT, consToMis, block, prune, revgap, min_length);  // for individual element
                cout<<"element "<<to_string(c)<<", number of base pair="<<to_string(bppc.GG)<<endl;
                if(filter) {
                  if(verbose) cerr << "filter: "<< c <<endl;
                  continue;
                }

                int tot = 0;
                if(bppc.idblk_count==0){
                    double nblk=(double)(bppc.GG - 15)/block;
                    int nblk2=(bppc.GG - 15)/block;
                    if( (nblk- nblk2)< ((double) 1.0/3.0)){
                        tot=nblk2;
                    }else{
                        tot=nblk2+1;
                    }
                    //tot=ceil((double)(bppc.GG - 15)/block);  //first bp length is always 15: len={0,15}. If change 15, change in Gibbs as well.
                    //have to make sure 15<=block. So better change 15 to block
                }else{
                    double nblk = (double)(bppc.GG - bppc.idblk_count-15)/block;
                    int nblk2 = (bppc.GG - bppc.idblk_count-15)/block;
                    if((nblk- nblk2)< ((double) 1.0/3.0)){
                        tot = 1+nblk2;
                    }else{
                        tot = 2+nblk2;
                    }
                    //tot=1+ceil((double)(bppc.GG - bppc.idblk_count-15)/block); 
                }
                // int tot=ceil((double)(bppc.GG - 15)/block); 
                //cout<<"tot ="<<tot<<endl;

                //null model
                cout<<"start null model\n";
                for (int iter = 0; iter <= tot; iter++)
                {
                    // cout<<"start initMCMC"<<endl;
                    bppc.initMCMC(iter, tot, bpp, 0, prune, false);

                    // cout<<"start Gibbs"<<endl;
                    bppc.Gibbs(iter, tot, bpp, out_Z0, output_path, output_path2, 0, true, sample_hyper, lrate_prop, grate_prop, WL); // Gibbs run to get Z for each element
                    if (bppc.verbose || bppc.failure)
                    {
                        bppc.Output_sampling(iter, output_path2, bpp, 0);
                        bppc.Output_tree(iter, output_path2, bpp, 0);
                    }
                }
                bppc.Output_init(output_path,output_path2,bpp,out_Z0, out_tree0, bppc.verbose); //sort rates!!, posterior median of nrate and crate; posterior mean of Z

                // //res model
                cout<<"start restricted model\n";
                for (int iter = 0; iter <= tot; iter++)
                {
                    bppc.initMCMC(iter, tot, bpp, 2, prune, false);
                    // cout<<"start Gibbs"<<endl;
                    bppc.Gibbs(iter, tot, bpp, out_Z2, output_path, output_path2, 2, true, sample_hyper, lrate_prop, grate_prop, WL);
                    if (bppc.verbose || bppc.failure)
                    {
                        bppc.Output_sampling(iter, output_path2, bpp, 1);
                        bppc.Output_tree(iter, output_path2, bpp, 1);
                    }
                }
                bppc.Output_init(output_path,output_path2,bpp,out_Z1, out_tree1, bppc.verbose);

                // full model
                cout<<"start full model\n";
                for (int iter = 0; iter <= tot; iter++)
                {
                    bppc.initMCMC(iter, tot, bpp, 1, prune, false); // not constrain log_prob_back
                    bppc.Gibbs(iter, tot, bpp, out_Z2, output_path, output_path2, 1, true, sample_hyper, lrate_prop, grate_prop, WL);
                    if (bppc.verbose || bppc.failure)
                    {
                        bppc.Output_sampling(iter, output_path2, bpp, 2);
                        bppc.Output_tree(iter, output_path2, bpp, 2);
                    }
                }
                bppc.Output_init(output_path,output_path2,bpp,out_Z2, out_tree2, bppc.verbose);
                
                cout << c << "\t" << bpp.log_liks_WL[0][c] <<"\t" <<  bpp.log_liks_WL[2][c] <<"\t" <<  bpp.log_liks_WL[1][c] <<endl;
                cout<<"\t" << bpp.log_liks_Z[0][c] << "\t" << bpp.log_liks_Z[2][c]<<"\t" << bpp.log_liks_Z[1][c] << endl;

            }catch (exception& e){
              cout << c << " Standard exception: " << e.what() << endl;
            }
      }

      bpp.Output_init(profile,output_path, ids);
    }

    out_Z0.close();
    out_Z1.close();
    out_Z2.close();
    out_hyper.close();
    out_lik.close();
    out_tree0.close();
    out_tree1.close();
    out_tree2.close();

    cout << endl << endl << "time used:  " << (time(NULL)-start)/60 << " min." << endl << endl;
    return 0;
}
