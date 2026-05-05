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


#include "../PhyloAcc-common/newick.h"
#include "../PhyloAcc-common/utils.h"
#include "../PhyloAcc-common/bpp_data.h"
#include "../PhyloAcc-common/bpp_io.h"
#include "../PhyloAcc-common/bpp_mcmc.h"
#include "../PhyloAcc-common/bpp_tree.h"
#include "../PhyloAcc-common/bpp_likelihood.h"
#include "../PhyloAcc-common/bpp_init.h"
#include "../PhyloAcc-common/bpp_active_output.h"
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
    phyloacc::InitializeTreeArrays(tree, N, -1, children, parent, distances);

    
//    distances[children[N-1][1]] += distances[children[N-1][0]];
//    distances[children[N-1][0]] = 0;
//    moveroot = children[N-1][0];
    
    
    //distances[83] += distances[42];  // modify distance for root!
    //distances[42] =0 ;
    
    
}


// try to match the phylogenetic profile and tree
void BPP::MatchProfAndTree(PhyloProf & _prof, PhyloTree & _tree)
{
    phyloacc::MatchProfileToTree(_prof, _tree);
}

void BPP::InitMCMC(int _num_burn, int _num_mcmc, int _num_thin)
{
    // init parameters
    num_burn = _num_burn;
    
    num_mcmc = _num_mcmc;
    num_thin = _num_thin;
    
    
    last_time = time(NULL);
    
    
    
    phyloacc::InitializeCommonMCMCStorage(
        C, N, ratio0, ratio1, ind_lrate, ind_lrate2, ind_grate,
        Max_Z, cur_Z, log_liks_null, log_liks_Z, log_liks_sgl, log_liks_resZ,
        log_liks_curZ, log_liks_propZ, MH_ratio_gain, MH_ratio_loss,
        cur_crate, cur_nrate, cur_lrate, cur_lrate2, cur_grate);
    
}

void BPP::sample_proposal(int iter, double & lrate_prop, double & grate_prop, ofstream & output)
{
    phyloacc::SampleProposal(RNG, iter, ind_lrate, ind_grate, vlr, vgr, nprior_a, nprior_b,
                             cprior_a, cprior_b, indel, indel2, lrate_prop, grate_prop, output);
}



double BPP::log_lik(vector< vector<vec> > & lambda, double _indel, double _indel2, int start1, int end1, vector<unsigned int> & v, double p)
{
    return phyloacc::ComputeLogLikelihood(
        lambda, _indel, _indel2, start1, end1, v, p, subtree, S, children, distances,
        [](const mat& log_x, const vec& log_y) { return BPP::log_multi(log_x, log_y); },
        [](const vec& log_y) { return BPP::log_exp_sum(log_y); });
}

void BPP::sample_hyperparam(int iter, vector<int> & ids, ofstream & output)
{
    phyloacc::SampleHyperparameters(RNG, iter, ids, cur_nrate, cur_crate, cur_lrate, cur_lrate2, cur_grate,
                                    nprior_a, nprior_b, cprior_a, cprior_b, prior_l_a, prior_l_b,
                                    prior_l2_a, prior_l2_b, prior_g_a, prior_g_b, output);
}

void BPP::getUppertree(int root, vector<int>& child, set<int> & visited_init)
{
    phyloacc::CollectUpperTreeNodes(root, child, parent, visited_init);
}

void BPP::getSubtree(int root, vector<int> & visited_init)
{
    phyloacc::CollectSubtreeNodes(root, children, visited_init);
}

void BPP::getSubtree(int root, set<int>& child, vector<int> & visited_init)
{
    phyloacc::CollectSubtreeNodesUntilChildren(root, children, child, visited_init);
}

void BPP::Output_init(PhyloProf & prof, string output_path, vector<int> & ids){
    phyloacc::WriteSTElemLikelihoodFile(output_path, prof.element_names, ids, log_liks_null, log_liks_resZ, log_liks_sgl, log_liks_Z);
    phyloacc::WriteElemZFiles(output_path, N, nodes_names, ids, Max_Z);
}


void BPP::Output_init0(PhyloProf & prof, ofstream& out_lik, vector<int> & ids){
    phyloacc::WriteInitLikelihoods(prof, out_lik, ids, log_liks_sgl, log_liks_Z);
}
