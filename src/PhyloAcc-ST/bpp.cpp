//
//  bpp.cpp
//  PhyloAcc
//
//  Created by hzr on 3/8/16.
//  Copyright © 2016 hzr. All rights reserved.
//

#include "bpp.hpp"
#include <armadillo>
#include <sys/types.h>
#include <dirent.h>
#include<queue>

#include <cmath>
#include <cassert>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ctype.h>


#include "newick.h"
#include "utils.h"
#include "bpp_c.hpp"


using namespace std;
using namespace arma;

// load the phylogenetic tree
void BPP::InitPhyloTree(PhyloTree & tree) //, double indel_pi), double indel, double indel2
{
   
//    cx_mat bvec;
//    cx_vec aval;
//   
//    if(num_base <= 4)  // no indel
//    {
//       
//    }else{
//        tree.subs_rate *= (1-indel);
//        mat B = ones<mat>(4,1) * indel;
//        mat C = ones<mat>(1,5) * indel2;
//        //tree.subs_rate.diag() -= indel;
//        tree.subs_rate.insert_cols(4, B);
//        tree.subs_rate.insert_rows(4, C);
//        //tree.subs_rate(4,4) = -0.04;
//        colvec c = sum(tree.subs_rate,1);
//        tree.subs_rate.diag() -=c;
//        
////        tree.pi.insert_rows(4,1);
////        tree.pi.head(4) *= 1 - indel_pi;
////        tree.pi[4] = indel_pi;
//        
//    }
//    
//    eig_gen(aval, bvec, tree.subs_rate);
//    eigenval = conv_to<mat>::from(aval);
//    eigenvec = conv_to<mat>::from(bvec).t();
//    eigeninv = inv(eigenvec);
//    submat = tree.subs_rate;
//    
//    //cout <<eigenval;
//    //cout <<"eigenvec: " << eigenvec;
//    //cout <<tree.subs_rate;
//    
//    mat a = null(tree.subs_rate.t());
//    
//    
//    pi = a/accu(a); //tree.pi;
//    //cout <<"pi: " <<  pi.t();
//    
//    log_pi = log(pi);
    
    submat = tree.subs_rate;
    children    = new int[N][2];
    parent      = new int[N];
    distances   = new double[N];
  
    
    for(int s=0; s<N; s++)
    {
        distances[s] = tree.distances[s];
        
    }
    
    
    
    for(int i=0; i<N; i++)
    {
        children[i][0] = -1;
        children[i][1] = -1;
        parent[i] = -1;
    }
    
    for(int i=0; i<N; i++)
    {
        int p = -1;
        for(int j=0; j<N; j++)
        {
            if (tree.dag[i][j])
            {
                p++;
                children[i][p] = j;
                parent[j] = i;
            }
        }
    }
    
    
//    distances[children[N-1][1]] += distances[children[N-1][0]];
//    distances[children[N-1][0]] = 0;
//    moveroot = children[N-1][0];
    
    
    //distances[83] += distances[42];  // modify distance for root!
    //distances[42] =0 ;
    
    
}


// try to match the phylogenetic profile and tree
void BPP::MatchProfAndTree(PhyloProf & _prof, PhyloTree & _tree)
{
    
    // try to match the species
    bool success_match = true;
    int S = _tree.S;
    int S2 = _prof.S;
    vector<int> reorder(S);  //each species in the tree where is in prof
    for(int s1=0; s1<S; s1++)
    {
        bool has_same_species = false;
        string sname1 = _tree.species_names[s1];
        for(int s2=0; s2<S2; s2++)
        {
            string sname2 = _prof.species_names[s2];
            //            cout << sname1 << " ? " << sname2 << endl;
            if (sname1 == sname2)
            {
                has_same_species = true;
                reorder[s1] = s2;
                break;
            }
        }
        if (!has_same_species)
        {
            cout << "No matrix species " << _prof.species_names[s1] << " found in tree." << endl;
            success_match = false;
            break;
        }
    }
    
    if (!success_match) // if cannot match literally
    {
        cout << endl << "The species in phylogenetic profile and tree cannot be matched literally:" << endl;
        cout << "The program will use the default mapping in data:" << endl;
        for(int s=0; s<S; s++)
            cout << "(" << _prof.species_names[s] << "\t=  " << _tree.species_names[s] << ")" << endl;
        cout << endl;
    }
    else                // if successully matched
    {
        cout << "The species in profile and tree match perfectly. Reorder the species in profile matrix by the tree." << endl << endl;
        vector<string> old_X = _prof.X;
        for(int s=0; s<S; s++)
        {
            int reorder_s = reorder[s];
            _prof.X[s] = old_X[reorder_s];
        }
        _prof.species_names = _tree.species_names;
    }
    
}

void BPP::InitMCMC(int _num_burn, int _num_mcmc, int _num_thin)
{
    // init parameters
    num_burn = _num_burn;
    
    num_mcmc = _num_mcmc;
    num_thin = _num_thin;
    
    
    last_time = time(NULL);
    
    
    
    // init MCMC sampling storage
    Max_Z = vector<vector < vector <int > >> (3,vector < vector <int >>  (C,vector <int > (N,0)));
    cur_Z = vector<vector < vector <int > >> (3,vector < vector <int >>  (C,vector <int > (N,0)));
    
    log_liks_null = vector <double>(C,0);
    log_liks_Z = vector<vector <double>>(3,vector <double>(C,0));
    log_liks_sgl = vector <double>(C,0);
    log_liks_resZ = vector <double>(C,0);
    
    log_liks_curZ = vector <double>(C,0);
    log_liks_propZ = vector <double>(C,0);
    MH_ratio_gain = vector <double>(C,0);
    MH_ratio_loss = vector <double>(C,0);
    
    
    
    cur_crate = vector <double>(C,ratio0);
    cur_nrate = vector <double>(C,ratio1);
    
    cur_lrate = vector <double>(C,ind_lrate);
    cur_lrate2 = vector <double>(C,ind_lrate2);
    cur_grate = vector <double>(C, ind_grate);
    
    
    
}

void BPP::sample_proposal(int iter, double & lrate_prop, double & grate_prop, ofstream & output)
{
    //indel_prop =gsl_ran_gamma(RNG,20 ,indel/20);
    //indel2_prop =gsl_ran_gamma(RNG,10 ,indel2/10);
    lrate_prop =gsl_ran_beta(RNG, ind_lrate * vlr, (1 - ind_lrate) *vlr); // let vlr == vgr!!
    grate_prop =gsl_ran_beta(RNG, ind_grate * vgr, (1 - ind_grate - ind_lrate) *vgr); //gsl_ran_gamma(RNG, vgr, ind_grate/vgr);
    grate_prop = grate_prop * (1 - lrate_prop);
   
//    submat.submat(0,0,3,3) *= (1 - indel_prop)/(1 - indel);
//    submat.col(4) = ones<mat>(5,1) * indel_prop;
//    submat.row(4) = ones<mat>(1,5) * indel2_prop;
//    colvec c = sum(submat,1);
//    submat.diag() -=c;
//    
//    cx_mat bvec;
//    cx_vec aval;
//    eig_gen(aval, bvec, submat);
//    eigenvalprop = conv_to<mat>::from(aval);
//    eigenvecprop = conv_to<mat>::from(bvec).t();
//    eigeninvprop = inv(eigenvecprop);
//    
//    mat a = null(submat.t());
//    piprop = a/accu(a);
//    //log_piprop = log(a/accu(a));
//    
//    //cout <<"proposed subs matrix: " << submat << endl;
//    //cout <<"piprop: " <<  piprop.t();
    
    output << iter << "\t"<< nprior_a<< "\t"<< nprior_b <<"\t"<< cprior_a << "\t"<< cprior_b << "\t"<< indel << "\t"<< indel2 << "\t"<<ind_grate<< "\t"<< ind_lrate <<endl;
}



double BPP::log_lik(vector< vector<vec> > & lambda, double _indel, double _indel2, int start1, int end1, vector<unsigned int> & v, double p)
{
    // compute loglik
    double result =0;
    mat x(2,2);
    
    int rr = *subtree.rbegin();
    // 1. sending the lambda msg from leaves bottom up through the network
    for(vector<int>::iterator it = subtree.begin(); it!=subtree.end(); it++) //int s=S; s<N; s++)
    {
        int s = *it;
        if(s<S) continue;
        int* p = children[s];
        for(int it = start1; it < end1; it++)  lambda[v[it]][s].fill(0);
       
        
        for(int cc=0;cc<2;cc++)
        {
            int chi = p[cc];
            assert(chi != -1);
            if(distances[chi]>0 )
            {
                double tt = (1 - exp(-(_indel + _indel2) * distances[chi]))/(_indel + _indel2);
                x.at(1,0) = _indel * tt;
                x.at(0,0) = 1 - x.at(1,0);
                
                x.at(0,1) = _indel2 * tt;
                x.at(1,1) = 1 - x.at(0,1);
                
                //cout << x;
                x = log(x);
                

            }
            else{
                x.fill(-INFINITY); //83, root
                x.diag().fill(0);
            }
            
            #pragma omp parallel for schedule (guided)
            for(int it = start1; it < end1; it++)  lambda[v[it]][s] +=  BPP::log_multi(x,lambda[v[it]][chi]);
        }
        
    }
    
    // 2. processing the distribution of root species
    
    for(int it = start1; it < end1; it++)
    {
        
//        if(lambda[v[it]][N-1][1] < -1e3)
//        {
//            cout <<v[it]<<": "<< lambda[v[it]][N-1].t();
//            cout <<children[N-1][0]<<": " << lambda[v[it]][children[N-1][0]].t();
//            cout <<children[N-1][1]<<": " << lambda[v[it]][children[N-1][1]].t();
//            
//        }
        lambda[v[it]][rr][0] += log(1-p); //N-1
        lambda[v[it]][rr][1] += log(p) ;
        result += BPP::log_exp_sum(lambda[v[it]][rr]);
    }
    
    return(result);
}

void BPP::sample_hyperparam(int iter, vector<int> & ids, ofstream & output) // recompute log_TM, double indel_prop, double indel2_prop,
{
// 0702 no sample rates
    //indepent MH to sample hyperparam of rates
    double p=1,r = 1; //hyperparam for shape
    double q=0.1,s = 0.1; // hyperparam for scale

    double vna = 100, vnb = 100, vca = 100, vcb = 100;
    double nprior_a_prop =gsl_ran_gamma(RNG, vna, nprior_a/vna);
    double cprior_a_prop =gsl_ran_gamma(RNG, vca, cprior_a/vca);

    double nprior_b_prop =gsl_ran_gamma(RNG, vnb, nprior_b/vnb);
    double cprior_b_prop =gsl_ran_gamma(RNG, vcb, cprior_b/vcb);

    //MH proposal
    double sum_r = 0;
    double log_prod_r = 0;
    //double var_r = 0;
    //for(vector<double>::iterator it = cur_nrate.begin(); it< cur_nrate.end(); it++)
    for(std::size_t i = 0; i < ids.size(); i++ )  // 0702
    {
        int c = ids[i];
        sum_r += cur_nrate[c];
        //var_r += pow(*it, 2);
        log_prod_r += log(cur_nrate[c]);
    }


    double M_ratio = (nprior_a_prop - 1) * (log(p) + log_prod_r) - (q + sum_r)/nprior_b_prop - (ids.size() + r)*lgamma(nprior_a_prop) - log(nprior_b_prop) * nprior_a_prop * (s + ids.size());
    M_ratio -= (nprior_a - 1) * (log(p) + log_prod_r) - (q + sum_r)/nprior_b - (ids.size() + r)*lgamma(nprior_a) - log(nprior_b) * nprior_a * (s + ids.size());

    double H_ratio = log(gsl_ran_gamma_pdf(nprior_a,vna,nprior_a_prop/vna)) - log(gsl_ran_gamma_pdf(nprior_a_prop,vna,nprior_a/vna)) + log(gsl_ran_gamma_pdf(nprior_b,vnb,nprior_b_prop/vnb)) - log(gsl_ran_gamma_pdf(nprior_b_prop,vnb,nprior_b/vnb));

    cout << "nrate_MH_ratio: " << M_ratio <<", " << H_ratio << ", " << nprior_a << ", " << nprior_a_prop << ", " << nprior_b << ", " << nprior_b_prop << endl;

    if(log(gsl_rng_uniform(RNG)) < M_ratio + H_ratio)
    {
        nprior_a = nprior_a_prop;
        nprior_b = nprior_b_prop;
    }


    sum_r = 0;
    log_prod_r = 0;
    //for(vector<double>::iterator it = cur_crate.begin(); it< cur_crate.end(); it++)
    for(std::size_t i = 0; i < ids.size(); i++ )  // 0702
    {
        int c = ids[i];
        sum_r += cur_crate[c];
        //var_r += pow(*it, 2);
        log_prod_r += log(cur_crate[c]);
    }

    M_ratio = (cprior_a_prop - 1) * (log(p) + log_prod_r) - (q + sum_r)/cprior_b_prop - (ids.size() + r)*lgamma(cprior_a_prop) - log(cprior_b_prop) * cprior_a_prop * (s + ids.size()) - ((cprior_a - 1) * (log(p) + log_prod_r) - (q + sum_r)/cprior_b - (ids.size() + r)*lgamma(cprior_a) - log(cprior_b) * cprior_a * (s + ids.size()));

    H_ratio = log(gsl_ran_gamma_pdf(cprior_a,vca,cprior_a_prop/vca)) - log(gsl_ran_gamma_pdf(cprior_a_prop,vca,cprior_a/vca)) + log(gsl_ran_gamma_pdf(cprior_b,vcb,cprior_b_prop/vcb)) - log(gsl_ran_gamma_pdf(cprior_b_prop,vcb,cprior_b/vcb));

    cout << "crate_MH_ratio: " << M_ratio + H_ratio << ", " << cprior_a << ", " << cprior_a_prop <<", " << cprior_b << ", " << cprior_b_prop << endl;

    if(log(gsl_rng_uniform(RNG)) < M_ratio + H_ratio)
    {
        cprior_a = cprior_a_prop;
        cprior_b = cprior_b_prop;
    }


  // sample hyperparameters of lrate and grate, exponential prior for prior_l_a and prior_l_b, prior_g_a and prior_g_b
    double u = gsl_rng_uniform(RNG)*(1.3 - 0.7) + 0.7; // generate uniform from (1/prop_n, prop_n);
    double prior_l_a_prop = prior_l_a *  u;
    
    u = gsl_rng_uniform(RNG)*(1.3 - 0.7) + 0.7; // generate uniform from (1/prop_n, prop_n);
    double prior_l_b_prop = prior_l_b *  u;
    double log_p = 0, log_pc = 0;
    for(std::size_t i = 0; i < ids.size(); i++ )  // 0702
    {
        int c = ids[i];
        log_p += log(cur_lrate[c]);
        log_pc += log(1 - cur_lrate[c]);
    }
    
    
    M_ratio = (prior_l_a_prop - prior_l_a) * (log_p - 1) + (prior_l_b_prop - prior_l_b) * (log_pc - 1);
    M_ratio += ids.size() * (gsl_sf_lnbeta(prior_l_a, prior_l_b) - gsl_sf_lnbeta(prior_l_a_prop, prior_l_b_prop));
    
    H_ratio = log(prior_l_a) - log(prior_l_a_prop) + log(prior_l_b) - log(prior_l_b_prop);
    
    cout << "lrate_MH_ratio: " << M_ratio + H_ratio << ", " << prior_l_a << ", " << prior_l_a_prop <<", " << prior_l_b << ", " << prior_l_b_prop  << endl;
    
    if(log(gsl_rng_uniform(RNG)) < M_ratio + H_ratio)
    {
        prior_l_a = prior_l_a_prop;
        prior_l_b = prior_l_b_prop;
    }
    
    // sample grate
    u = gsl_rng_uniform(RNG)*(1.3 - 0.7) + 0.7; // generate uniform from (1/prop_n, prop_n);
    double prior_g_a_prop = prior_g_a *  u;
    
    u = gsl_rng_uniform(RNG)*(1.3 - 0.7) + 0.7; // generate uniform from (1/prop_n, prop_n);
    double prior_g_b_prop = prior_g_b *  u;
    log_p = 0; log_pc = 0;
    for(std::size_t i = 0; i < ids.size(); i++ )  // 0702
    {
        int c = ids[i];
        log_p += log(cur_grate[c]);
        log_pc += log(1 - cur_grate[c]);
    }
    
    
    M_ratio = (prior_g_a_prop - prior_g_a) * (log_p - 1) + (prior_g_b_prop - prior_g_b) * (log_pc - 1);
    M_ratio += ids.size() * (gsl_sf_lnbeta(prior_g_a, prior_g_b) - gsl_sf_lnbeta(prior_g_a_prop, prior_g_b_prop));
    
    H_ratio = log(prior_g_a) - log(prior_g_a_prop) + log(prior_g_b) - log(prior_g_b_prop);
    
    cout << "grate_MH_ratio: " << M_ratio + H_ratio << ", " << prior_g_a << ", " << prior_g_a_prop <<", " << prior_g_b << ", " << prior_g_b_prop  << endl;
    
    if(log(gsl_rng_uniform(RNG)) < M_ratio + H_ratio)
    {
        prior_g_a = prior_g_a_prop;
        prior_g_b = prior_g_b_prop;
    }
    
     output << iter << "\t"<< nprior_a<< "\t"<< nprior_b <<"\t"<< cprior_a << "\t"<< cprior_b << "\t"<< prior_l_a << "\t"<< prior_l_b << "\t"<< prior_g_a << "\t"<< prior_g_b <<endl;
    
   
}








void BPP::getUppertree(int root, vector<int>& child, set<int> & visited_init)  // include root!
{
    for(vector<int>::iterator it = child.begin(); it!=child.end(); it++)
    {
        int p = *it;
       while(p!=root)
       {
           visited_init.insert(p);
           p = parent[p];
       }
        
        
    }
    
    visited_init.insert(root);
    
    
    
}





void BPP::getSubtree(int root, vector<int> & visited_init)  // traverse from root to children, include root
{
    
    
    int j = root;
    
    //cout << nodes_names[j]<<"\t";
    
    if(children[j][0]!=-1)
    {
        
        for(int chi=0;chi<2;chi++)
        {
            getSubtree(children[j][chi], visited_init);
            
            
        }
        
    }
   
    visited_init.push_back(j);
    
}

void BPP::getSubtree(int root, set<int>& child, vector<int> & visited_init)  // traverse from root, stop at children, 74 & 64; do include 1-S!
{
    
    
    int j = root;
    
    //cout << nodes_names[j]<<"\t";
    if(child.find(j) != child.end())
    {
        
        visited_init.push_back(j);
        
        return;
    }
    
    if(children[j][0]!=-1)
    {
        
        for(int chi=0;chi<2;chi++)
        {
            getSubtree(children[j][chi], child, visited_init);
            
            
        }
        
      //  visited_init.push_back(j);
        
    }
    
     visited_init.push_back(j);
       
    
}

void BPP::Output_init(PhyloProf & prof, string output_path, vector<int> & ids){
    
        string outpath_elem = output_path+ "_elem_lik.txt";
        ofstream out_lik(outpath_elem.c_str());
        out_lik.precision(8);
        out_lik << "No.\tID\tloglik_Null\tloglik_Acc\tloglik_Full\tlogBF1\tlogBF2\tloglik_Max_M0\tloglik_Max_M1\tloglik_Max_M2"<<endl;
        //for(int cc=0; cc<C;cc++)
        for(vector<int>::iterator it = ids.begin(); it !=ids.end(); it++)
        {
            int cc = *it;
            out_lik <<cc << "\t" << prof.element_names[cc] << "\t" << log_liks_null[cc] <<"\t"  <<log_liks_resZ[cc] <<"\t"  <<log_liks_sgl[cc]<< "\t";
            out_lik << log_liks_resZ[cc] -  log_liks_null[cc] << "\t" << log_liks_resZ[cc] -  log_liks_sgl[cc];
            //for(int r=0;r<3;r++) 
            out_lik <<"\t" <<log_liks_Z[0][cc] << "\t" <<log_liks_Z[2][cc]<<"\t" <<log_liks_Z[1][cc];
            out_lik << endl;
        }
        
        out_lik.close();
        
    ofstream out_z;
    for(int r =0;r<3;r++)
    {
            if(r == 2)
            {
                outpath_elem = output_path+"_M" +to_string(1) + "_elem_Z.txt";
            }else if(r == 1)
            {
                outpath_elem = output_path+"_M" +to_string(2) + "_elem_Z.txt";
            }else{
                outpath_elem = output_path+"_M" +to_string(0) + "_elem_Z.txt";
            }
            
            out_z.open(outpath_elem.c_str());
            out_z<<"No."; 
            for(int s =0 ;s<N;s++){  // header: species name
                out_z<< "\t" << nodes_names[s] ;
            }
            out_z <<endl;
            
            for(vector<int>::iterator it = ids.begin(); it !=ids.end(); it++)
            {
                int c = *it;
                out_z<<c; 
                for(int s=0; s<N;s++)
                    out_z <<"\t" << Max_Z[r][c][s];
                out_z <<endl;
            }
        out_z.close();
    }
    
    
}


void BPP::Output_init0(PhyloProf & prof, ofstream& out_lik, vector<int> & ids){
    
    //for(int cc=0; cc<C;cc++)
    for(vector<int>::iterator it = ids.begin(); it !=ids.end(); it++)
    {
        int cc = *it;
        out_lik <<cc << "\t" << prof.element_names[cc] << "\t"  <<log_liks_sgl[cc]<< "\t"<<log_liks_Z[1][cc];
        out_lik << endl;
    }
    
    //out_lik.close();
    
    
}




