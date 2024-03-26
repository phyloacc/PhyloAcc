//
//  main.cpp
//  PhyloAcc
//
//  Created by hzr on 3/8/16.
//  Copyright © 2016 hzr. All rights reserved.
//

/////////////////////////////////////////////////////////////////

#include <dirent.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iomanip>
//#include <omp.h>
#include <string> 
#include <armadillo>
#include "profile.h"
#include "newick.h"
#include "bpp.hpp"
#include "bpp_c.hpp"
#include "utils.h"
#include <math.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace arma;

/////////////////////////////////////////////////////////////////

// parameters, input and output files paths
string params_path;
string phytree_path;
string align_path;
string output_path="";
string output_path2="";
string segment_path;
string id_path="";
string result_prefix="test";

//string refspecies;
string outgroup;
string targetspecies;
string conservegroup; // can't be missing in more than 50%
double conserve_prop = 0.8;

int num_thread = 1;

// running parameters
int num_burn = 200;         // num of burn-in updates
int num_mcmc = 800;         // num of MCMC updates,10
int num_thin = 500;    // num of updates between two samples, adaptive frequency
int num_chain; // outer loop of updates Q matrix and hyperparameter of substitution rates


double prep_lrate = 0.3;
double prep_lrate2 = 0.1;
double prep_grate = 0.5; // initalization

double prior_lrate2_a = 1,prior_lrate2_b = 9 ; // beta prior for lrate2, 0.5
double prior_lrate_a = 1 ,prior_lrate_b = 9 ; // beta prior for lrate, 2
double prior_grate_a = 3,prior_grate_b = 1 ; // beta prior for grate


double ratio0 = 0.5; //initial conserved rate
double ratio1 = 1; // initial accelerated rate
double missing_thres = 0.8;

double nprior_a = 10, nprior_b = 0.2;  //around 1
double cprior_a = 5, cprior_b = 0.04;  //around ratio
int ropt = 1;
double cub = 1;
double nlb = 0.6;

int batch = -1 ;
int seed = 5;
double indel = 0;
double indel2;
bool sample_indel = 0;
bool sample_hyper;
char gapchar = '-';
bool verbose = 0;
double consToMis = 0.01;
bool prune=0;
double revgap=1;
int min_length = 50;

/////////////////////////////////////////////////////////////////
// Functions

bool stringToBool(const std::string& s) {
    return s == "true" || s == "1";
}

// load the program parameters
void LoadParams(int argc, char* argv[])
{
    cout << "Loading input data and running parameters......" << endl;
    
    if (argc > 1)
    {
        params_path = string(argv[1]);
    }
    else
    {
        //params_path = "/Users/zhu/PhyloAcc/code/PhyloAcc-test-data/phyloacc-test-config.txt";  //params2
        cout << "Please specify the path of the parameter file." << endl;
        exit(1);
    }
    // Get the name of the parameters file, or error out if none is provided

    cout << "Loading program configurations from " << params_path << "......" <<endl;
    
    const int BUFF_SIZE = 4096; // max line length
    char line_buff[BUFF_SIZE];

    vector<string> param_line;
    string word;
    int num_words;
    int line_num = 1;
    // Variables for reading the parameter file

    ifstream in_params(params_path.c_str());
    if (!in_params)
    {
        cerr << "Cannot open the parameters file: " << params_path.c_str() << endl;
        exit(1);
    }
    // Check if the parameter file exists and can be opened

    cout << "Reading parameters......" << endl;
    while(!in_params.eof())
    {
        param_line.clear();
        in_params.getline(line_buff, BUFF_SIZE);
        istringstream line_stream(line_buff);
        //string tmp; line_stream >> tmp;
        
        while(line_stream >> word) {
            param_line.push_back(word);
        }
        // Read the line into line_buff, then parse it into param_line

        num_words = param_line.size();
        // Count the number of words in the line

        // for(const auto &word : param_line) {
        //     cout << word << ' ';
        // }
        // cout << '\n';
        // cout << num_words << endl;
        // Debugging

        if (param_line[0][0] == '#' || num_words == 0)
            continue;
        // Skip comment lines and empty lines in the config file

        if (num_words != 2) {
            cerr << endl << "Line " << line_num << " in the parameter file is not formatted correctly: " << line_buff << endl;
            cerr << "Each line should contain a parameter and a value separated by a space and no other whitespace." << endl << endl;
            exit(1);
        }
        // Check that the line is formatted correctly

        // input and output file paths
        if (param_line[0]=="PHYTREE_FILE")
            phytree_path = param_line[1];
        else if (param_line[0]=="ALIGN_FILE")
            align_path = param_line[1];
        else if (param_line[0]=="SEG_FILE")
            segment_path = param_line[1];
        else if (param_line[0]=="ID_FILE")
            id_path = param_line[1];
        else if (param_line[0]=="BATCH")
            batch = stoi(param_line[1]);
        else if (param_line[0]=="RESULT_FOLDER")
            output_path = param_line[1];
        else if (param_line[0]=="PREFIX")
            result_prefix = param_line[1];
                
        //else if (tmp=="RESULT_INDIV")
        //    line_stream >> output_path2;
        
        else if (param_line[0]=="SEED")
            seed = stoi(param_line[1]);
        else if (param_line[0]=="INIT_CONSERVE_RATE")
            ratio0 = stod(param_line[1]);
        else if (param_line[0]=="INIT_ACCE_RATE")
            ratio1 = stod(param_line[1]);
        else if (param_line[0]=="CONSERVE_PRIOR_A")
            cprior_a = stod(param_line[1]);
        else if (param_line[0]=="CONSERVE_PRIOR_B")
            cprior_b = stod(param_line[1]);
        else if (param_line[0]=="ACCE_PRIOR_A")
            nprior_a = stod(param_line[1]);
        else if (param_line[0]=="ACCE_PRIOR_B")
            nprior_b = stod(param_line[1]);
        else if (param_line[0]=="ROPT")
            ropt = stoi(param_line[1]);
        else if (param_line[0]=="CUB")
            cub = stod(param_line[1]);
        else if (param_line[0]=="NLB")
            nlb = stod(param_line[1]);

        // running parameters
        else if (param_line[0]=="BURNIN")
            num_burn = stoi(param_line[1]);
        else if (param_line[0]=="MCMC")
            num_mcmc = stoi(param_line[1]);
        else if (param_line[0]=="ADAPT_FREQ")
            num_thin = stoi(param_line[1]);
        else if (param_line[0]=="INIT_LRATE")
            prep_lrate = stod(param_line[1]);
        else if (param_line[0]=="INIT_LRATE2")
            prep_lrate2 = stod(param_line[1]);
        else if (param_line[0]=="INIT_GRATE")
            prep_grate = stod(param_line[1]);
        else if (param_line[0]=="HYPER_LRATE_A")
            prior_lrate_a = stod(param_line[1]);
        else if (param_line[0]=="HYPER_LRATE_B")
            prior_lrate_b = stod(param_line[1]);
        else if (param_line[0]=="HYPER_GRATE_A")
            prior_grate_a = stod(param_line[1]);
        else if (param_line[0]=="HYPER_GRATE_B")
            prior_grate_b = stod(param_line[1]);
        else if (param_line[0]=="HYPER_LRATE2_A")
            prior_lrate2_a = stod(param_line[1]); // if it's zero, assume dollo
        else if (param_line[0]=="HYPER_LRATE2_B")
            prior_lrate2_b = stod(param_line[1]);
        else if (param_line[0]=="CHAIN")
            num_chain = stoi(param_line[1]);
        
        // constraint
        else if (param_line[0] == "OUTGROUP")
            outgroup = param_line[1];
        else if (param_line[0] == "TARGETSPECIES")
            targetspecies = param_line[1];
        else if (param_line[0] == "CONSERVE")
            conservegroup = param_line[1];
        else if (param_line[0] == "CONSERVE_PROP")
            conserve_prop = stod(param_line[1]);
        else if (param_line[0] == "GAP_PROP")
            missing_thres = stod(param_line[1]);
        //else if (tmp == "REF")
         //   line_stream >> refspecies;
        else if (param_line[0] == "CONSTOMIS")
            consToMis = stod(param_line[1]);
        
        
        // treat indel as additional character
        else if (param_line[0] == "GAPCHAR")
            gapchar = param_line[1][0];
        else if (param_line[0] == "PRUNE_TREE")
            prune = stringToBool(param_line[1]);
        else if (param_line[0] == "TRIM_GAP_PERCENT")
            revgap = stod(param_line[1]);
        else if (param_line[0] == "MIN_LEN")
            min_length = stoi(param_line[1]);
        else if (param_line[0] == "INDEL") // not used
            indel = stoi(param_line[1]);
        else if (param_line[0] == "INDEL2") // not used
            indel2 = stoi(param_line[1]);
        else if(param_line[0] == "SAMPLE_INDEL")  // not used
            sample_indel = stringToBool(param_line[1]);
        else if(param_line[0] == "SAMPLE_HYPER")
            sample_hyper = stringToBool(param_line[1]);
        else if(param_line[0] == "VERBOSE")
            verbose = stringToBool(param_line[1]);
        else if(param_line[0] == "NUM_THREAD")
            num_thread = stoi(param_line[1]);
        else
            cout << "Unknown parameter: " << param_line[0] << ", skipping..." << endl;

        line_num++;
        // Increment the line number
    }
    
    // trimming file names
    phytree_path = strutils::trim(phytree_path, " \"\t\n");
    align_path = strutils::trim(align_path, " \"\t\n");
    output_path  = strutils::trim(output_path,  " \"\t\n");
    //output_path2  = strutils::trim(output_path2,  " \"\t\n");
    segment_path = strutils::trim(segment_path, " \"\t\n");
    
    if(prior_lrate2_a==0) prep_lrate2 = 0;
}

////////////////////

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

////////////////////

void DispParams(PhyloProf profile, int seed)
{
    double mean_seg_size = 0;
    for(unsigned int c=0; c<profile.C; c++)
        mean_seg_size += (double)(profile.element_pos[c][1] - profile.element_pos[c][0]) / profile.C;
    cout << "# total length = " << profile.G << " (" << profile.C << ")" << ". # Species = " << profile.S << ". # elements = " << profile.C << ". Mean gene set size = " << mean_seg_size << "." << endl;
    cout << "# Burn-ins = " << num_burn << ". # MCMC Updates = " << num_mcmc << ". # adaptive frequency = " << num_thin << ".  RND SEED = " << seed << "." << endl ; //
    cout << "# Threads = " << num_thread << endl << endl;
}

/////////////////////////////////////////////////////////////////
// Main

int main(int argc, char* argv[])
{
    time_t start = time(NULL);
    
    cout << std::fixed << setprecision(4);
    
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
    
    // load the phylogenetic tree
    PhyloTree phytree = LoadPhyloTree(phytree_path);
    
    // init and display the running parameters
    //InitParams(profile.G);
    DispParams(profile, seed);
    
    // create and init the BPP object
    //int pC = 500;  // only read in some elements for testing
    BPP bpp(0, profile, phytree, output_path, targetspecies, outgroup, conserve_prop, conservegroup, ratio0, ratio1, ropt, cub, nlb, nprior_a, nprior_b, cprior_a, cprior_b, seed, prep_grate, prep_lrate, prep_lrate2, prior_grate_a, prior_grate_b,prior_lrate_a, prior_lrate_b,prior_lrate2_a, prior_lrate2_b,  indel, indel2, missing_thres, sample_indel);  //c=1 test run first element
    
    // remove profile?
    //profile.~PhyloProf();
    
    //initialize the MCMC sampling
    bpp.InitMCMC(num_burn, num_mcmc, num_thin);
    
    output_path = output_path + "/" + result_prefix ;
    output_path2 = output_path;
    string outpath_Z0 = output_path + "_rate_postZ_M" +to_string(0) +".txt";  //null
    string outpath_Z1 = output_path + "_rate_postZ_M" +to_string(2) +".txt";  //full
    string outpath_Z2 = output_path + "_rate_postZ_M" +to_string(1) +".txt";  // M1
    string outpath_hyper = output_path+"_hyper.txt";
    // Output file names    
    
    ofstream out_hyper(outpath_hyper.c_str());
    out_hyper << "iter\tnprior_a\tnprior_b\tcprior_a\tcprior_b\tprior_l_a\tprior_l_b\tprior_g_a\tprior_g_b\n";
    out_hyper << 0 << "\t"<< nprior_a<< "\t"<< nprior_b <<"\t"<< cprior_a << "\t"<< cprior_b << "\t"<< prior_lrate_a << "\t"<< prior_lrate_b << "\t"<< prior_grate_a << "\t"<< prior_grate_b <<endl;
    // Hyperparameter output headers

    ofstream out_lik;
    if(sample_hyper) {
        string outpath_elem = output_path+ "_elem_lik.txt";
        out_lik.open(outpath_elem.c_str());
        out_lik.precision(8);
        out_lik << "No.\tID\tloglik_all\tloglik_Max"<<endl;
    }
    // Elem likelihood output headers
    
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
    // Z file headers

    double lrate_prop = 0.5, grate_prop = 0.5;
    
    vector<int> ids;
    if(id_path=="")
    {
        if(batch==-1)
        {
          for(int c =0;c<bpp.C;c++)
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

    cout << ids.size() << " of elements to be computed" << endl; 
    //gsl_set_error_handler_off();
    for(int iter =0; iter<num_chain; iter++)
    {
        cout << "Running MCMC chain " << iter +1 << " ..." << endl;
        // generate next indel parameter, nprior_a, nprior_b, etc...
        //bpp.sample_proposal(iter, lrate_prop, grate_prop,out_hyper);
       
        
        // Gibbs sampling
        #pragma omp parallel for schedule (guided) num_threads(num_thread)
        for(std::size_t i = 0; i < ids.size(); i++ ) 
        {
            int c = ids[i];
            bool filter = false;
          
            // null model
            try{
            BPP_C bppc(c, profile, bpp, gapchar, missing_thres, filter, verbose, consToMis, prune, revgap, min_length);  // for individual element
            if(filter) {
              if(verbose) cerr << "filter: "<< c <<endl;
              continue;
            }
            
            if(!sample_hyper)
            {
                bppc.initMCMC(iter,bpp,0);  //not constrain log_prob_back
                bppc.Gibbs(iter,bpp,out_Z0,output_path,output_path2,0,true,sample_hyper, lrate_prop, grate_prop);  // Gibbs run to get Z for each element, lrate_prop & grate_prop not used 
                //if(!bppc.failure)
                //{
                bppc.Eval2(bpp,0);
                if(bppc.verbose || bppc.failure) bppc.Output_sampling(iter, output_path2, bpp, 0);
                bppc.Output_init(output_path,output_path2,bpp,out_Z0, 0); //sort rates!!, posterior median of nrate and crate; posterior mean of Z
               // }
                // res model, crate by null model
                bppc.initMCMC(iter,bpp,2);  //not constrain log_prob_back
                bppc.Gibbs(iter,bpp,out_Z2,output_path,output_path2,2,true,sample_hyper, lrate_prop, grate_prop);  // Gibbs run to get Z for each element
                //if(!bppc.failure)
                //{
                bppc.Eval2(bpp,2);
                if(bppc.verbose || bppc.failure) bppc.Output_sampling(iter, output_path2, bpp, 1);
                bppc.Output_init(output_path,output_path2,bpp,out_Z2, 2); //sort rates!!
                //}
            }
            
            // full model, nrate, crate by res model
            bppc.initMCMC(iter,bpp,1);  //not constrain log_prob_back
            bppc.Gibbs(iter, bpp,out_Z1,output_path,output_path2,1, true, sample_hyper, lrate_prop, grate_prop);  // Gibbs run to get Z for each element
            //if(!bppc.failure)
            //{
            bppc.Eval2(bpp,1);
            if(bppc.verbose || bppc.failure) bppc.Output_sampling(iter, output_path2, bpp, 2);
            bppc.Output_init(output_path,output_path2,bpp,out_Z1, 1); //sort rates!!
            //}

            }catch (exception& e){
              cout << c << " Standard exception: " << e.what() << endl;
            }
        }
        
        // sample indel, nprior_a, nprior_b, etc...
        try{
          if(sample_hyper) {
            bpp.sample_hyperparam(iter, ids, out_hyper);
            bpp.Output_init0(profile,out_lik, ids);
          }else{
            bpp.Output_init(profile,output_path, ids);
          }
        }catch (exception& e){
              cout << " Standard exception: " << e.what() << endl;
        }
    }

    
        out_Z0.close();
        out_Z1.close();
        out_Z2.close();
        out_hyper.close();
        out_lik.close();
    
    
    
    cout << endl << endl << "time used:  " << (time(NULL)-start)/60 << " min." << endl << endl;
    
    
    return 0;
}

/////////////////////////////////////////////////////////////////