//
//  bpp_c.cpp
//  PhyloAcc_init2
//
//  Created by hzr on 4/19/16.
//  Copyright © 2016 hzr. All rights reserved.
//

#include "bpp_c.hpp"
#include "bpp.hpp"
#include <gsl/gsl_errno.h>

//struct Cmp
//{
//    bool operator ()(const pair<int, double> &a, const pair<int, double> &b)
//    {
//        return a.second <= b.second;
//    }
//};

void BPP_C::getSubtree(int root, vector<int> & visited_init)  // traverse from root to children, include root
{
    int j = root;

    for(int chi=0;chi<2;chi++)
    {
        if(children2[j][chi]!=-1)
            getSubtree(children2[j][chi], visited_init);
    }
    visited_init.push_back(j);
}


void BPP_C::getUppertree_missing(BPP & bpp)  // set children2 to -1 if missing & in child
{
    vector<bool> temp_missing = vector<bool>(N, false);
    
   for(set<int>::iterator it= bpp.upper.begin(); it!=bpp.upper.end(); it++) //lineages leading to outgroups 
   {
       int s = *it;
       if(s < S && missing[s]) temp_missing[s] = true;
   }
    
    // reset missing only both children in upper & missing
    for(int s=S; s<N; s++)
    {
        int* p = bpp.children[s];
        if(temp_missing[p[0]] && temp_missing[p[1]]) {
            temp_missing[s] = true;
        }
    }
    
    for(set<int>::iterator it= bpp.upper.begin(); it!=bpp.upper.end(); it++)
    {
        int s = *it;
        if(s < S) continue;
        int* p = bpp.children[s];
        for(int cc =0 ;cc < 2; cc++)
        {
            if(temp_missing[p[cc]])
            {
                children2[s][cc] = -1;
            }
        }
    }    
}


void BPP_C::simulate(BPP& bpp, PhyloProf & profile, char gapchar, bool prune)
{
    Z = vector<int>(N,1);
    bpp.cur_Z[0][0]=vector<int> (N,1);
    
    //get common ancestors of target_species
    vector<bool> visited = vector<bool>(N);
    vector<int> nontarget = vector<int>();
    nontarget.reserve(bpp.conservedgroup.size() + bpp.outgroup.size());
    nontarget.insert(nontarget.end(), bpp.outgroup.begin(), bpp.outgroup.end());
    nontarget.insert(nontarget.end(), bpp.conservedgroup.begin(), bpp.conservedgroup.end());
    
    getUpdateNode(nontarget, visited, bpp); // visited = false for all upstream of nontarget_species
    //Han*: whole lineage leading to a nontarget, visited=False. For rest: visited=True
    
    for(int s = 0; s < N; s++){
        if(visited[s]){
            Z[s] = 2;
            bpp.cur_Z[0][0][s]=2; //Han* for output Z
        }  
    }
    
    double nrate=gsl_ran_gamma(RNG,bpp.nprior_a,bpp.nprior_b);
    double crate = gsl_ran_gamma(RNG, bpp.cprior_a, bpp.cprior_b);
    if(crate>nrate){
        nrate=crate*1.2;
    }

    bpp.cur_nrate[CC] = nrate;
    bpp.cur_crate[CC] = crate;

    pi[0]=(double) gsl_ran_beta(RNG,prior_dir_param[0],prior_dir_param[1])/2;
    pi[3]=pi[0];
    pi[1]=(double) (1-2*pi[0])/2;
    pi[2]=pi[1];
    bpp.cur_pi[0][CC][0]=pi[0]; //For writing output 
    log_pi=log(pi); 
    
    getQmat(pi,inst_rate);
    
    // sample gene tree
    if(prune)
    {
        set<int> alls;
        for(int s=0; s<N; s++) alls.insert(s);
        gtree->initTree(missing, alls, bpp);
    }else{
        gtree->initTree(missing, bpp.upper, bpp);
    }
    
    for(int s= S; s<N; s++)
    {
        for(map<int, int>::iterator it = gtree->parent_gene2[s].begin(); it!= gtree->parent_gene2[s].end(); it++ )
        {
            gtree->Tg[s][it->first] = vector<int>(GG, -2);

        }
        
    }
   
    gtree->InitTg(GG, bpp, nodes);
    gtree->Simulate_Tg(GG, bpp, nodes, Z, bpp.cur_nrate[CC], bpp.cur_crate[CC], log_pi, eigenvec, eigenval, eigeninv);

    int st = bpp.element_start[CC];
    for(int s =0; s<S; s++)
    {
        for(int g = 0; g < GG; g++){
            switch (gtree->Tg[s][s][g])
            {
                case 0:
                    profile.X[s][g+st] = 'a'; break;
                case 1:
                    profile.X[s][g+st] = 'c'; break;
                case 2:
                    profile.X[s][g+st] = 'g'; break;
                case 3:
                    profile.X[s][g+st] = 't'; break;
                case 4:
                    profile.X[s][g+st] = gapchar; break;
                default:
                    profile.X[s][g+st] = 'n';
             
            }
        }
    }
    std::stringstream buffer;
    gtree->printTree(gtree->root, bpp, buffer);
    bpp.genetrees[0][CC] = buffer.str();
}

//assume GTR model
void BPP_C::getQmat(vec & statR, vec & instR)
{
    submat(0, 0) = -(instR(0) * statR(1) + instR(1) * statR(2) + instR(2) * statR(3));
    submat(0, 1) = instR(0) * statR(1);
    submat(0, 2) = instR(1) * statR(2);
    submat(0, 3) = instR(2) * statR(3);
    submat(1, 0) = instR(0) * statR(0);
    submat(1, 1) = -(instR(0) * statR(0) + instR(3) * statR(2) + instR(4) * statR(3));
    submat(1, 2) = instR(3) * statR(2);
    submat(1, 3) = instR(4) * statR(3);
    submat(2, 0) = instR(1) * statR(0);
    submat(2, 1) = instR(3) * statR(1);
    submat(2, 2) = -(instR(1) * statR(0) + instR(3) * statR(1) + instR(5) * statR(3));
    submat(2, 3) = instR(5) * statR(3);
    submat(3, 0) = instR(2) * statR(0);
    submat(3, 1) = instR(4) * statR(1);
    submat(3, 2) = instR(5) * statR(2);
    submat(3, 3) = -(instR(2) * statR(0) + instR(4) * statR(1) + instR(5) * statR(2));
    
    cx_mat bvec;
    cx_vec aval;
    eig_gen(aval, bvec, submat);
    eigenval=conv_to<mat>::from(aval);
    eigenvec=conv_to<mat>::from(bvec).t();
    eigeninv=inv(eigenvec);
}

void BPP_C::initMCMC(int iter, int max_iter, BPP&bpp, int resZ, bool prune, bool fixtree)
{

    failure = false;
    prop_c = 1.2;  //(double) 1/GG * 5; //
    prop_n = (double) 1/GG * 50;  // 20, probably change to larger for smaller sample size at beginning
    Z = bpp.cur_Z[resZ][CC];
    trace_genetree.clear();
    accept_n_rate = 0;
    accept_c_rate = 0;
    
    if(iter == 0)
    {       
        //init a new tree
        if(!fixtree)
        {
            if(prune)
            {
                set<int> alls;
                for(int s=0; s<N; s++) alls.insert(s);
                gtree->initTree_Sptop(missing, alls,bpp);
            }else{
                gtree->initTree_Sptop(missing, bpp.upper,bpp);
                //gtree->printSptree(bpp,1);
            }
        }else{
            gtree->initTree(bpp.element_tree[CC], bpp);
        }
        
        trace_n_rate[0] = bpp.cur_nrate[CC]; //1;
        trace_c_rate[0] = bpp.cur_crate[CC]; //ratio;

        trace_l_rate[0] = bpp.cur_lrate[CC];
        trace_l2_rate[0] = bpp.cur_lrate2[CC];
        trace_g_rate[0] = bpp.cur_grate[CC];

        std::fill(trace_GTtopChg.begin(),trace_GTtopChg.end(),0);
        
    }else{
        int tmp_m=num_burn+num_mcmc-1;

        trace_n_rate[0] = trace_n_rate[tmp_m]; //1;
        trace_c_rate[0] = trace_c_rate[tmp_m]; //ratio;

        trace_l_rate[0] = trace_l_rate[tmp_m];
        trace_l2_rate[0] = trace_l2_rate[tmp_m];
        trace_g_rate[0] = trace_g_rate[tmp_m];

        std::fill(trace_GTtopChg.begin(),trace_GTtopChg.end(),0);
        trace_pi[0]=trace_pi[tmp_m]; 
        if(prune)
        {
            set<int> alls;
            for(int s=0; s<N; s++) alls.insert(s);
            gtree->initTree_Sptop(missing, alls, bpp);
        }else{
            gtree->initTree_Sptop(missing, bpp.upper, bpp);
        }
   
    }
    
    for(int s=S; s<N; s++) 
    {
        gtree->lambda[s].clear();
        gtree->Tg[s].clear();
    }
    

    fixZ = vector<int>(N,-1);
    for(vector<int>:: iterator it = upper_c.begin(); it < upper_c.end(); it++)
    {
        fixZ[*it] = 1;
    }

    fixZ[root] =1;  // originally: 0... (0: 0, 1: 0/1, -1: all)

    MaxLoglik = -INFINITY;

    log_prob_back = vector<vec> (N, zeros<vec>(3));

    //with prior, null model, all branches except outgroup are conserved, only have to sample Z  (0/1) of outgroup
    if(resZ==2)
    {
        for(vector<int>::iterator it = upper_conserve_c.begin(); it < upper_conserve_c.end(); it++)
        {
            int s= *it;
            fixZ[s] = 1;
        }
    }
    else if(resZ==0)
    {
        for(vector<int>::iterator it = nodes.begin();it<nodes.end()-1;it++)  //don't include root, subtree doesn't have the root edge
        {
            int s = *it;
            fixZ[s] = 1;
        }
   }

    for(vector<int>::iterator it = nodes.begin(); it <nodes.end(); it++)
    {
        int s = *it;
        if(fixZ[s] == 1)
        {
            log_TM_Int[s](0,1) = log(0);  // set once, never changed
            log_TM_Int[s](1,1) = 0;
            log_TM_Int[s](2,1) = log(0);

            log_TM_Int[s](0,0) = log(1 - trace_g_rate[0]);
            log_TM_Int[s](1,0) = log(trace_g_rate[0]); //grate=alpha in PHI. prob of bkgrd -> conserv
            log_TM_Int[s](2,0) = log(0);
            
            // won't use log_TM_Int(,2)
        
        }else{
            log_TM_Int[s](0,0) = log(1 - trace_g_rate[0]);
            log_TM_Int[s](1,0) = log(trace_g_rate[0]);
            log_TM_Int[s](2,0) = log(0) ;
            
            double y = 1 - trace_l_rate[0];
            log_TM_Int[s](0,1) = log(0); // set once, never changed
            log_TM_Int[s](1,1) = log(y);
            log_TM_Int[s](2,1) = log(1-y);
            
            log_TM_Int[s](0,2) = log(0);
            log_TM_Int[s](1,2) = log(trace_l2_rate[0]);
            log_TM_Int[s](2,2) = log(1 - trace_l2_rate[0]);
        }
    }
    //printSptree(bpp);
}


// ofstream & outZ, no use!,
// iter start from 0 to blocks = 10, return the ratio
void BPP_C::Gibbs(int iter, int max_iter, BPP &bpp, ofstream & outZ, string output_path,string output_path2,int resZ, bool UpR, bool UpHyper, double lrate_prop, double grate_prop, bool WL)
{

    vector<int> changedZ;
    vector<bool> visited = vector<bool>(N,false);
    vec logp_Z = zeros<vec>(2);
    double accept_tree = 0;
    double accept_tree_bn=0; //acceptance rate in burnin.
    //double accept_tree_bn2=0;
    double accept_tree_br = 0, total_sample_tree = 0;
    //double bothAccept=0; //acceptance rate of branchlen of the topo changed part.
    double accept11 = 0, total_br = 0; //acceptance rate of branchlen over the tree
    double MHdelta=1; //for branchLen update over the tree
    double MHdelta_pi=0.1; //unif proposal half width for pi
    double accept_pi =0, total_pi = 0;

    int indictor = 0; // start from more complex model
    vec ind_prop = zeros<vec>(2);
    vec logNormalConst = zeros<vec>(2), log_lik_old = zeros<vec>(2);
    vec momentum = zeros<vec>(2);
    double learning_rate = 1;
    double beta = 0.9;
    vector<int> lens = vector<int>(2, 0);
    bool recordtree=false;

    //get prob of sampling potential gene tree branch
    int tot  = gtree->var_br_node.size(); 
    int ngnode = gtree->gene_nodes.size();
    vec prob_br_sample(tot);    
    if(tot<(ngnode-1)){ //root node is not in var_br_node
        prob_br_sample.fill((double) 1.0/tot); //if only uncertain for a few branches, then give equal prob of sampling these branches.
    }else{ //if all br can be sampled, then sample reversely proportional to sp pa length.
        double sum_prob=0;
        for(int tmp_i=0; tmp_i<tot; tmp_i++){
            prob_br_sample[tmp_i]=gtree->prob_var_node[tmp_i];
            sum_prob +=prob_br_sample[tmp_i];
        }
        prob_br_sample /=sum_prob;
    }

    if(WL)
    {   
        if(iter==0 && idblk_count>0){
            lens[0]=idblk_count;
            lens[1]=0;
        }else{ // if(iter != max_iter)
            if(idblk_count>0){
                lens[0]=idblk_count+(iter)* GG_block;
                if(lens[0] > GG || iter==max_iter) lens[0] = GG;
                lens[1] = idblk_count+(iter-1)* GG_block; //have to make sure GG_block>=15;
                if(lens[1] <= 0) lens[1] = 0; 
            }else{
                lens[0]=(iter)* GG_block+15;
                if(lens[0] > GG || iter==max_iter) lens[0] = GG;
                lens[1] = (iter-1)* GG_block+15; //have to make sure GG_block>=15;
                if(lens[1] <= 0) lens[1] = 0;
            }
        } 
    }else{
        lens[0] = GG;
    }
   
   if(lens[0]==GG) recordtree=true;
    
    double logNorm = 0;

    //initial lambda and missing patterns of Tg
    gtree->GG = lens[0];
    for(int s= S; s<N; s++)
    {
        for(map<int, int>::iterator it = gtree->parent_gene2[s].begin(); it!= gtree->parent_gene2[s].end(); it++ )
        {          
            gtree->lambda[s][it->first].resize(gtree->GG, zeros<mat>(bpp.num_base,3));
            gtree->Tg[s][it->first].resize(gtree->GG, -2); //resize(n, val): n-newsize, val-initialized value
        }
        
    }
    vector<double> obsACGT=vector<double>(4,0.0);
    for(int g=0; g<lens[0]; g++){
        for(int s=0; s<S; s++){
            int a = gtree->Tg[s][s][g];
            if(a>=0 && a<=3) obsACGT[a]+=1;
        }
    }

    //gtree->InitTg(lens[0], bpp, nodes); //11Nov need this? everytime initMCMC, gtree is re-initiliazed, Tg are cleared.
    gtree->Update_Tg(lens[0], visited, bpp, false, nodes, Z, trace_n_rate[0], trace_c_rate[0], log_pi, eigenvec, eigenval, eigeninv);

    GTree gt_init(N, S);
    gtree->copyto(lens[0], gt_init); // lambda filled by initial values, Tg is -2/4

    double lpz = gtree->get_logpZ1(bpp, Z, consToMis, nconsToMis); //P(missing | Z, GT), changed from get_logpZ 0517

    if(iter !=0)
    {
        logp_Z[0] = logp_Z[1] = lpz; //+priorZ;
    }else{
        logp_Z[1] =0;// priorZ;
        logp_Z[0] = lpz ; //+priorZ;
    }
    
    double logP_tree = gtree->priorTree(bpp);
    double priorZ = priorP_Z(bpp);
    
    int new_burn=num_burn;
    int new_mcmc=num_mcmc;

    //cout<<"Gibbs iter="<<iter<<endl;
    for(m=0; m<(new_burn+new_mcmc); m++)
    {        
        //cout<<"m="<<m<<endl;
        if(m==0)
        {
            trace_loglik[m] = 0;
            // calculate marginal likelihood,
            for(int g=0;g<lens[0];g++)
            {
                double temp = BPP::log_exp_sum(gtree->lambda[N-1][gtree->root][g].col(2) + log_pi); //proportional to P(Yg | rest)
                trace_loglik[m] += temp; 
                if(g <lens[1])
                    log_lik_old[1] += temp; 
            }
            log_lik_old[0] =trace_loglik[m];  

            // initialize logNormalConst
            for(int ii =0; ii<2;ii++)
                logNormalConst[ii] = log_lik_old[ii] + logp_Z[ii]; //-priorZ; //P(Y,Missing |Z, rest param),
        }else{
            trace_loglik[m] = log_lik_old[0];
        }
     
        if(WL)
        {
            vec prob = BPP::log_sample_norm(log_lik_old + logp_Z - logNormalConst); //Han: return normalized prob
            unsigned int n[2];
            gsl_ran_multinomial(RNG, 2, 1, prob.memptr(),n);
            int nn;
            for(nn=0; nn<2;nn++)
            {
                if(n[nn]>0) break;
            }
            indictor = nn; // WL which dist to sample.
            if(m<(new_burn+new_mcmc-1)) trace_indicator[m+1]=indictor;
          
            momentum = beta*momentum - learning_rate*prob;

            ind_prop[indictor] ++;
            logNormalConst -= momentum;
           
            double normalize_result = BPP::log_exp_sum(logNormalConst);
            logNormalConst -= normalize_result;

            //adjust the learning rate if necessary
            double rr = ind_prop[0] / ind_prop[1];
            if(rr > 2.0/3.0 && rr <3.0/2.0){
                learning_rate = 1/(1 + 1/learning_rate);
                //cout<<"Model res="<<resZ<<", iter"<<iter<<", m="<<m<<". rr="<<rr<<"; eta="<<learning_rate<<endl;
                ind_prop.fill(0);
            }

            if(m == new_burn+new_mcmc - 1)
            {
                double logNorm1 =0;
                logNorm = 0;
                for(std::size_t i =new_burn; i< new_burn+new_mcmc;i++)
                {
                    logNorm += trace_logNormratio[i];
                    if(i >= new_burn + ceil(new_mcmc/2)) logNorm1 += trace_logNormratio[i];
                }
                logNorm /=new_mcmc;
                logNorm1 /= new_mcmc-ceil(new_mcmc/2);
                //cout<<"Model res="<<resZ<<", iter"<<iter<<", m="<<m<<". logNorm1:="<<logNorm1<<"; 2:="<<logNorm<<endl;
                
                if(abs(logNorm1 - logNorm) > 0.5) { //*maxr.second - *maxr.first
                    new_burn = trace_logNormratio.size() - new_mcmc;
                }
            }
            
        }
        
        trace_logNormratio[m] = logNormalConst[0] - logNormalConst[1];

        MonitorChain(m, iter, max_iter, bpp, ind_prop[0] / ind_prop[1], logp_Z[0] + priorZ + logP_tree , resZ, recordtree);  //+ prior_Z

        if(failure || (m == new_burn+new_mcmc-1))  break;

        // regenerate gt from species tree if accept rate is too low
        if(idblk_count==0 || iter>0){
        if(m<= new_burn && total_sample_tree == 100){ 
            if(accept_tree_bn/total_sample_tree < 0.05){
                gtree->copyfrom(lens[0],gt_init,bpp.num_base); // initial Tg == -2 here!
                gtree->Update_Tg(lens[0], visited, bpp, false, nodes, Z, trace_n_rate[0], trace_c_rate[0], log_pi, eigenvec, eigenval, eigeninv);
                double lpz = gtree->get_logpZ1(bpp,Z, consToMis, nconsToMis);
                if(iter !=0)
                {
                    logp_Z[0] = logp_Z[1] = lpz; 
                }else{
                    logp_Z[1] =0;
                    logp_Z[0] = lpz ; 
                }
                for (int g = 0; g < lens[0]; g++)
                {
                    double temp = BPP::log_exp_sum(gtree->lambda[N - 1][gtree->root][g].col(2) + log_pi);
                    if (g < lens[1])
                        log_lik_old[1] += temp;
                    log_lik_old[0] += temp;
                }
                
                //cout << "regenerate gt from species, accept rate: " << accept_tree_bn/total_sample_tree;
                //cout << "log_lik_old: " << log_lik_old[0] << ", " << log_lik_old[1] << "; lpz, missing pattern: " << lpz << endl;
            }
           
            accept_tree_bn = 0;
            total_sample_tree = 0;
        }else if(m == new_burn){
            accept_tree = 0;
            total_sample_tree = 0;
            accept_tree_bn = 0;
        }
        }

        if((m + 1) % num_thin==0){   
            int res_pi=sample_pi(m,indictor,lens,resZ,pi[0],log_lik_old,MHdelta_pi,bpp, obsACGT,0);
            accept_pi += res_pi; 
            total_pi += 1;
            
            if (total_pi == (10*num_thin) && (m+1) < new_burn){ //if acceptRate too low, re-sample a new pi
                
                if(accept_pi/total_pi < 0.1) {
                    res_pi = sample_pi(m, indictor, lens, resZ, pi[0],log_lik_old,MHdelta_pi,bpp, obsACGT,1);
                }
                accept_pi = 0;
                total_pi = 0;
            }        
        }else{
            for(int i=0; i<4; i++){
                trace_pi[m+1][i]=pi[i];
            }
        }

        // sample H_m|Z_m,r_m
        if(idblk_count==0 || iter>0 || (m + 1) % num_thin==0){
        gtree->Update_Tg(lens[indictor], visited, bpp, true, nodes, Z, trace_n_rate[m], trace_c_rate[m], log_pi, eigenvec, eigenval, eigeninv);
        // sample Z_(m+1)|H_m,r_m,t_m
        changedZ = Update_Z(lens[indictor], bpp);

        getUpdateNode(changedZ,visited,bpp);  // reset visited to true
        for (int i = 0; i < N; i++)  visited[i] = false;
       
        if(changedZ.size()>0){
            lpz = gtree->get_logpZ1(bpp, Z, consToMis, nconsToMis);
            priorZ=priorP_Z(bpp);
            if (iter != 0)
            {
                logp_Z[0] = logp_Z[1] = lpz ; 
            }
            else
            {
                logp_Z[1] = 0; 
                logp_Z[0] = lpz; 
            }
            gtree->Update_Tg(lens[indictor], visited, bpp, false, nodes, Z, trace_n_rate[m], trace_c_rate[m], log_pi, eigenvec, eigenval, eigeninv);
        }

        log_lik_old[0] = log_lik_old[1] = 0;
        for (int g = 0; g < lens[indictor]; g++){
            double temp = BPP::log_exp_sum(gtree->lambda[N - 1][gtree->root][g].col(2) + log_pi);
            if (g < lens[1])
                log_lik_old[1] += temp;
            log_lik_old[0] += temp;
        }
        } // reduce sampling freq if in idblk
      
        // sample r_(m+1)|Z_(m+1), r_m
        if( UpR && (m + 1) % num_thin==0 )  // update transition rate
        {
            vector<bool> visited_neut = vector<bool>(N,true);
            vector<bool> visited_cons = vector<bool>(N,true);

            getUpdateNode(1,visited_neut, bpp);  //all neut nodes, no lambda
            getUpdateNode(0,visited_cons, bpp);  //all conserved nodes

            // sample neutral and conserved rate by MH, c has to be after n!
            //int M =1, bool adaptive = true, double adaptive_factor = 0.5
            trace_n_rate[m+1] = sample_rate(indictor,iter, max_iter,  lens, resZ, trace_n_rate[m],true, visited_neut, log_lik_old,bpp, 1, true, 0.5);  //partial lambda changed, given Z, integrate history
            trace_c_rate[m+1] = sample_rate(indictor,iter, max_iter,  lens, resZ, trace_c_rate[m],false, visited_cons, log_lik_old,bpp, 1, true, 0.5);
        }else{
            trace_n_rate[m+1] = trace_n_rate[m];
            trace_c_rate[m+1] = trace_c_rate[m];
        }
       
        if(idblk_count==0 || iter>0){ //if identical blks, keep gt=sp, no update.
        if(iter <= 5 || (m+1) % num_thin==0)
        {
            if(lens[indictor]>0){   
            bool tempaccept=false;
            GTree gt_org(N, S);  // copy a genetree
            gtree->copyto(lens[indictor], gt_org);

            unsigned int n[tot];
            int branchTop;
            bool nonRootSample = true;
            while (nonRootSample)
            {
                gsl_ran_multinomial(RNG, tot, 1, prob_br_sample.memptr(), n);
                int nn;
                for (nn = 0; nn < tot; nn++)
                {
                    if (n[nn] > 0)
                        break;
                }
                if (gt_org.var_br_node[nn] != gt_org.root)
                {
                    branchTop = gt_org.var_br_node[nn];
                    nonRootSample = false;
                }
            }

                //printSptree(bpp);
                tempaccept = gt_org.Sample_tree2(branchTop, indictor, bpp,Z, trace_n_rate[m+1], trace_c_rate[m+1], consToMis, nconsToMis, lens, log_lik_old, logp_Z, eigenvec, eigenval, eigeninv, log_pi);
            
            total_sample_tree++;
           
            if(tempaccept) {
                gtree->copyfrom(lens[indictor], gt_org, bpp.num_base); // clear lambda and Tg from lens[indictor] -> GG
                logP_tree = gtree->priorTree(bpp);
                
                if(m>=new_burn) accept_tree += 1; else accept_tree_bn += 1;
                trace_GTtopChg[m+1]=1;
            }

            double accept1=0.0;
            if(m%num_thin==0){
            for (int b = 0; b <gtree->gene_nodes.size() ; b++)
            {   
                if (gtree->gene_nodes[b] >= S){ //sample internal node's kid's branch
                    GTree gt_org(N, S);
                    gtree->copyto(lens[indictor], gt_org);
                    bool tempacceptbr = gt_org.Sample_BranchLen(MHdelta,gtree->gene_nodes[b], indictor, bpp, Z, trace_n_rate[m + 1], trace_c_rate[m + 1], consToMis, nconsToMis, lens, log_lik_old, logp_Z, eigenvec, eigenval, eigeninv, logP_tree, log_pi);
                    if(tempacceptbr){
                        gtree->copyfrom(lens[indictor], gt_org, bpp.num_base);
                        logP_tree = gtree->priorTree(bpp);
                        accept1 = accept1 + 1;
                        accept_tree_br += 1;
                    }
                }
            }
            if(m>=new_burn) {
                accept11=accept11+accept1/(gtree->gene_nodes.size()-1); //later to change denominator
                total_br += 1;
            }

            }//sample br length 
        } //no sample if len[indictor]=0. 
            
        } //step size: sample every numThin steps or m<5
        }//no sampling if in idblk

        if(idblk_count==0 || iter>0 || (m + 1) % num_thin==0){
        if (indictor == 1){
            gtree->Update_Tg(lens[0], visited, bpp, false, nodes, Z, trace_n_rate[m + 1], trace_c_rate[m + 1], log_pi, eigenvec, eigenval, eigeninv, lens[1]);
            for (int g = lens[1]; g < lens[0]; g++){
                log_lik_old[0] += BPP::log_exp_sum(gtree->lambda[N - 1][gtree->root][g].col(2) + log_pi);
            }
        }
        } //reduce sampling freq in idblk;


        //cout<<"print trees after all param update"<<endl;
        //gtree -> printSptree(bpp,1);

        if(idblk_count==0 || iter>0 || (m + 1) % num_thin==0){
            sample_transition(trace_g_rate[m + 1] ,trace_l_rate[m + 1] ,trace_l2_rate[m + 1], bpp );
        }else{
            trace_g_rate[m+1]=trace_g_rate[m];
            trace_l_rate[m+1]=trace_l_rate[m];
            trace_l2_rate[m+1]=trace_l2_rate[m];
        }

        if(total_br == 100){ 
            if(accept11/total_br < 0.1){  // 3.5?
                MHdelta=max(MHdelta/2,0.1);
                //cout<<"delta="<<MHdelta<<endl;
            }else if(accept11/total_br > 0.8){
                MHdelta=min(MHdelta*2,5.0);
                //cout<<"delta="<<MHdelta<<endl;
            }
            //cout<<"m = " <<m << " accept rate of sampling branch len: "<<accept11/total_br<<"\tdelta = "<<MHdelta<<endl; //        if((m+1)%1000==0 && m>=new_burn)
            accept11 = 0;
            total_br = 0;
        }

        if (m == 0){ 
            bpp.log_mle[resZ][CC] = log_lik_old[0] + logp_Z[0];
        }else if (m >= new_burn & (log_lik_old[0] + logp_Z[0]) > bpp.log_mle[resZ][CC]){
            bpp.log_mle[resZ][CC] = log_lik_old[0] + logp_Z[0];
        }
    }
    
    bpp.log_liks_Z[resZ][CC] = MaxLoglik;
    bpp.Max_Z[resZ][CC] = Max_Z;
    bpp.genetrees[resZ][CC] = Max_GT; 
    bpp.cur_Z[resZ][CC] = Z; 
    //Han*: save running max
    bpp.cur_pi[resZ][CC]=Max_pi; 


    if(UpHyper) {        
        bpp.cur_crate[CC] = trace_c_rate[m];
        bpp.cur_nrate[CC] = trace_n_rate[m]; 
        
        bpp.cur_grate[CC] = trace_g_rate[m];
        bpp.cur_lrate[CC] = trace_l_rate[m];
        bpp.cur_lrate2[CC] = trace_l2_rate[m];
    }

    
     #pragma omp critical
    {
        if(verbose) {
            cout <<CC <<": " << Max_m << ", " << MaxLoglik << ", "<< trace_c_rate[Max_m] << ", " << trace_n_rate[Max_m];
            cout << ", "<< trace_g_rate[Max_m] << ", " << trace_l_rate[Max_m] << ", " << trace_l2_rate[Max_m] << endl;
            for(int s=0;s<N;s++ )
            {
                if(Max_Z[s]!=1 && !missing[s])
                {
                    cout << s << ": " <<Max_Z[s] << ", ";
                }
            }
            cout << endl << "gap: ";

            for(int s=0;s<N;s++ )
            {
                if(missing[s])
                {
                    cout << s << ": " <<Max_Z[s] << ", ";
                }
            }
            cout << endl << endl;
        }
    }
    
    if(WL) bpp.log_liks_WL[resZ][CC] += logNorm;

    //pass Mx est to next M_x+1 if possible
    if(resZ==0){
        bpp.cur_pi[2][CC]=bpp.cur_pi[0][CC];
    }else if(resZ==2){
        bpp.cur_Z[1][CC]=bpp.cur_Z[2][CC];
        bpp.cur_pi[1][CC]=bpp.cur_pi[2][CC];
    }
    //cout<<"iter"<<iter<<": overall AR of sample tree after burnin = "<<(double) accept_tree/total_sample_tree<<". totBrChg="<<accept_tree_br/(m+1)* num_thin/(gtree->gene_nodes.size()-S)<<endl;  
}

double BPP_C::priorP_Z(BPP & bpp)
{
    // compute P(Z|C)
    mat nZ = zeros(3,3);
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
    {
        int p = bpp.parent[*it];
        //if(Z[p] < 2)
        //{
            nZ(Z[*it],Z[p]) += 1;
        //}
    }
    double priorZ = gsl_sf_lnbeta(prior_g_a + nZ(1,0) + nZ(2,0), prior_g_b + nZ(0,0)) - gsl_sf_lnbeta(prior_g_a, prior_g_b);

    for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
    {
        if(fixZ[*it] == 1)
        {
            int p = bpp.parent[*it];
            assert(p!=N);
            //if(Z[p] < 2)
            //{
                nZ(Z[*it],Z[p]) -= 1;
            //}
        }
    }

    priorZ += gsl_sf_lnbeta(prior_l_a + nZ(2,1), prior_l_b + nZ(1,1)) - gsl_sf_lnbeta(prior_l_a, prior_l_b);
    if(prior_l2_a > 0) priorZ += gsl_sf_lnbeta(prior_l2_a + nZ(1,2), prior_l2_b + nZ(2,2)) - gsl_sf_lnbeta(prior_l2_a, prior_l2_b);

    //root
    priorZ += prior_z.at(Z[N-1]);
    
    return(priorZ);
}


void BPP_C::sample_transition( double  & gr, double  & lr, double  & lr2, BPP & bpp)
{

    mat nZ = zeros(3,3); // record number of Z transitions
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
    {
        int p = bpp.parent[*it];
        //if(Z[p] < 2)
        //{
            nZ(Z[*it],Z[p]) += 1;
        //}


    }

    //cout << "Count matrix: " << nZ << endl;
    gr = gsl_ran_beta(RNG, prior_g_a + nZ(1,0) + nZ(2,0) , prior_g_b + nZ(0,0));

    //for(vector<int>:: iterator it = upper_c.begin(); it < upper_c.end() -1;it++)
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
    {
        if(fixZ[*it] == 1)
        {
            int p = bpp.parent[*it];
            assert(p!=N);
            //if(Z[p] < 2)
            //{
                nZ(Z[*it],Z[p]) -= 1;
            //}
        }
    }

    lr = gsl_ran_beta(RNG, prior_l_a + nZ(2,1), prior_l_b + nZ(1,1));
    if(prior_l2_a == 0)
    {
        lr2 = 0;
    }else{
        lr2 = gsl_ran_beta(RNG, prior_l2_a + nZ(1,2), prior_l2_b + nZ(2,2));
    }

    for(vector<int>::iterator it = nodes.begin(); it <nodes.end(); it++)
    {
        int s = *it;

        if(fixZ[*it] == 1)
        {
            log_TM_Int[*it](1,1) = 0;
            log_TM_Int[*it](2,1) = log(0);

            log_TM_Int[*it](0,0) = log(1 - gr);
            log_TM_Int[*it](1,0) = log(gr);
            log_TM_Int[*it](2,0) = log(0);

        }else{
            log_TM_Int[s](0,0) = log(1 - gr);
            log_TM_Int[s](1,0) = log(gr);

            double y = 1 - lr;
            log_TM_Int[s](1,1) = log(y);
            log_TM_Int[s](2,1) = log(1-y);
            
            log_TM_Int[s](2,2) = log(1 - lr2);
            log_TM_Int[s](1,2) = log(lr2);
        }

    }

}

// only update probability of nodes above changedZ
void BPP_C::getUpdateNode(vector<int> changedZ, vector<bool> & visited_init, BPP & bpp) //changedZ from small to large (bottom to top)
{
    for(int i =0;i<N;i++)
        visited_init[i] = true;

    if(changedZ.size()==0) return;
    for(vector<int>::iterator it = changedZ.end()-1 ; it >= changedZ.begin(); --it){
        int j = *it;
        //j = bpp.parent[j];

        while(j!=N)
        {
            if(!visited_init[j]) break;
            visited_init[j] = false;
            //for(int g=0; g<GG; g++) lambda[g][j].zeros();  // only fathers of changedZ!

            j = bpp.parent[j];
            assert(j!=-1);

        }
    }

}

void BPP_C::MonitorChain(int  m,int iter, int max_iter, BPP &bpp, double ind_prop, const double add_loglik,const int resZ, bool recordtree)  //calculate P(X|Z, TM(r) ),
//void BPP_C::MonitorChain(int  m, BPP &bpp, double ind_prop, const double add_loglik,const int resZ, bool recordtree, int lensC)  //Han*: Debug
{


    for(int s=0; s<N;s++)
    {
        trace_Z[m][s] = Z[s];
    }

   // P(X, M,  r, T|Z) = P(X, M |r, Z, T ) * P(r)  * P(T)
    // trace_loglik P(Y|r, Z, T)
   // add_loglik = P(T) * P(M | Z, T )
   trace_full_loglik[m] = trace_loglik[m] + add_loglik ; //P(Y,Missing, Z, G| hyper etc)

   trace_full_loglik[m] += log(gsl_ran_gamma_pdf(trace_c_rate[m],bpp.cprior_a,bpp.cprior_b)); //P(Y,Miss,Z,G,rc |rest)

    if(resZ !=0)
   {    //P(Y,Miss,Z,G,rc, rn |rest)
        trace_full_loglik[m] += log(gsl_ran_gamma_pdf(trace_n_rate[m],bpp.nprior_a,bpp.nprior_b));

    }

    //P(Y,Miss,Z,G,rc, rn, pi |rest)
    trace_full_loglik[m]+=log(gsl_ran_beta_pdf(2*trace_pi[m][0],prior_dir_param[0],prior_dir_param[1]));
    
    //recordtree
    int new_burn=num_burn;

    if(m>=new_burn && recordtree)
    {
        int minss = -1;
        // save children_gene
        vector<vector<int>> temp_children=vector<vector<int>> (N, vector<int> (2,0));
        for(int s = 0; s< N; s++)
        {
            temp_children[s][0] = gtree->children_gene[s][0];
            temp_children[s][1] = gtree->children_gene[s][1];
        }
        
        gtree->getGeneNodes(gtree->root, minss); //minss=label of the left child, i.e., the child with the smaller label.
        vector<int> parent_relabel = vector<int>(N, -1);
        vector<double> branchlen = vector<double>(N, 0);
        int currentS = S;
        gtree->printTree2(gtree->root, currentS, parent_relabel, branchlen);
        for(int s = 0; s< N; s++)
        {
            gtree->children_gene[s][0] = temp_children[s][0];
            gtree->children_gene[s][1] = temp_children[s][1];
        }
        
        string tree_str ="";
        for(vector<int>::iterator it = parent_relabel.begin(); it!= parent_relabel.end(); it++)
        {
            tree_str += to_string(*it) + "," ;  // only record topology
        }
        map<string, simpletree>::iterator it = trace_genetree.find(tree_str);
        
        if(it!=trace_genetree.end())
        {
            it->second.count++;
            for(int s = 0; s< it->second.root+1; s++)
            {
                it->second.distances[s] += branchlen[s];
            }
        }else{
            int N_new=currentS +1 ;
            simpletree st = {1, vector<double>(N_new), vector<vector<int>>(N_new, vector<int>(2)), vector<int>(N_new), gtree->root, vector<string>(N_new)}; //Han*: add node names 
            for(int s = 0; s< N_new; s++)
            {
                st.parent_gene[s] = parent_relabel[s];
                st.distances[s] = branchlen[s];
                for(int cc =0; cc < 2; cc++) st.children_gene[s][cc] = -1;
                
            }
            for(int s = 0; s< N_new; s++)
            {
                int p = parent_relabel[s];
                if(p<N && p != -1)
                {
                    if(st.children_gene[p][0] == -1)
                    {
                        st.children_gene[p][0] = s;
                    }else{
                        st.children_gene[p][1] = s;
                    }
                }else if(p == N)
                {
                    st.root = s;
                }
            }
            //update node_names
            for(int s=0; s<S; s++){
                st.node_names[s]=to_string(s);
            }
            for(int s=S; s<N_new; s++){
                st.node_names[s]="("+st.node_names[st.children_gene[s][0]]+","+st.node_names[st.children_gene[s][1]]+")";
            }
            trace_genetree[tree_str] = st; 
        }
    }
    
    if(m>=new_burn && trace_full_loglik[m]>MaxLoglik)
    {
        MaxLoglik = trace_full_loglik[m];
      
        Max_Z = Z;
        Max_m = m;
        Max_pi=trace_pi[m];
            
        std::stringstream buffer;
        gtree->printTree(gtree->root, bpp, buffer);
        Max_GT = buffer.str();
        
    }

    if( m>=new_burn && (m%5000==0) && verbose){  // && verbose
        cout <<"iter"<<m<<", "<<CC <<": " << trace_loglik[m] <<", " << trace_full_loglik[m]<< ", " << trace_n_rate[m] << ", " << trace_c_rate[m] <<", " <<  prop_c <<", " << prop_n <<", ";
        cout << trace_g_rate[m] <<  ", " << trace_l_rate[m] << ", " << trace_l2_rate[m] << endl;
        cout << trace_logNormratio[m] <<", " << ind_prop << endl;
        if(recordtree) cout << "number of unique tree: " << trace_genetree.size() << endl;
    }
}


void BPP_C::getEmission(int len, BPP & bpp)
{
    
    if(len==0)
    {
        for(int s =0 ;s <N; s++)
        {
            log_emission[s][0] = 0;
            log_emission[s][1] = 0;
            log_emission[s][2] = 0;
        }
        return;
    }
    vector<double> rates; rates.push_back(1); rates.push_back(trace_c_rate[m]);rates.push_back(trace_n_rate[m]);
    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end();it++) //N-1
    {
        int ss = *it;
        int sp = bpp.parent[ss];
        double height_ps = sp < N? bpp.heights[sp] : INFINITY;
        
        vector<double> emissions = vector<double>(3, 0); //

        //go through all gtree nodes (exclude root), P(Missing /Not | Z), missing from the root side
        /*for(vector<int>::iterator it2 = gtree->temp_coal[ss].begin(); it2 != gtree->temp_coal[ss].end(); it2++){
            int coalNode=*it2;
            if(!gtree->missing_gene[coalNode]){ //assume missing at top, to be consistent with get_logPZ1.
                for(int cc=0; cc<2; cc++){
                    int cl=gtree->children_gene[coalNode][cc];
                    if(gtree->missing_gene[cl]){
                        emissions[0] += log(nconsToMis);  //logP(Missing_kids(gnode ins s) |Z_s ) for each Z=0:2 scenario.
                        emissions[1] += log(consToMis);
                        emissions[2] += log(nconsToMis);
                    }else{
                        emissions[0] += log(1 - nconsToMis);
                        emissions[1] += log(1 - consToMis);
                        emissions[2] += log(1 - nconsToMis);
                    }
                }
            }
        } */

        //missing from the bottom
        for(map<int,int>::iterator it2 = gtree->parent_gene2[ss].begin(); it2 != gtree->parent_gene2[ss].end(); it2++)
        {
            int cl = it2->first;
            if(cl<S){
                if (gtree->missing_gene[cl])
                {
                    emissions[0] += log(nconsToMis); //logP(Missing_kids(gnode ins s) |Z_s ) for each Z=0:2 scenario.
                    emissions[1] += log(consToMis);
                    emissions[2] += log(nconsToMis);
                }
                else
                {
                    emissions[0] += log(1 - nconsToMis);
                    emissions[1] += log(1 - consToMis);
                    emissions[2] += log(1 - nconsToMis);
                }
            }else{
                if(gtree->heights_gene[cl]>=bpp.heights[ss]){
                    if (gtree->missing_gene[cl])
                    {
                        emissions[0] += log(nconsToMis); //logP(Missing_kids(gnode ins s) |Z_s ) for each Z=0:2 scenario.
                        emissions[1] += log(consToMis);
                        emissions[2] += log(nconsToMis);
                    }
                    else
                    {
                        emissions[0] += log(1 - nconsToMis);
                        emissions[1] += log(1 - consToMis);
                        emissions[2] += log(1 - nconsToMis);
                    }
                }
            }
        } //emission: if Z=0,1,2, what are the P(miss/not |Z) prob.
    
        //sum over all g, integrate internal Tg
        vector<mat> logTMs = vector<mat>(2);
        map<int, vector<vec> > templambda;
        
        for(map<int, int>::iterator it = gtree->parent_gene2[ss].begin(); it!= gtree->parent_gene2[ss].end(); it++)
        {
            int gs = it->first; //kid
            if(gtree->heights_gene[gs] <= bpp.heights[ss])
            {
                templambda[gs] = vector<vec>(len, zeros<vec>(bpp.num_base));
                if(gtree->missing_gene[gs]) continue;

                for(int g=0;g < len;g++)
                {
                    if(gtree->Tg[ss].at(gs)[g] >= bpp.num_base) continue;
                    templambda[gs][g].fill(LOG_ZERO);
                    if(gtree->Tg[ss][gs][g] == -1){  //for r,w, etc
                       templambda[gs][g] = gtree->lambda[ss][gs][g].col(2);
                    }else if(ss >= S)
                    {
                        templambda[gs][g][gtree->Tg[ss][gs][g]] = 0;
                    }
                }
            }
        }

        //for leaf: below is not executed.So leaf nodes are not contributing to emission.
        for(int rr = 0; rr < 3; rr++)
        {
            for(vector<int>::iterator it2 = gtree->temp_coal[ss].begin(); it2 < gtree->temp_coal[ss].end(); it2++)
            {
                int gs = *it2;
                templambda[gs] = vector<vec>(len, zeros<vec>(bpp.num_base));
                vector<int> childs;
                for(int cc=0;cc<2;cc++)
                {
                    int chi = gtree->children_gene[gs][cc];
                    if(gtree->missing_gene[chi]) continue;  // for speed
                    childs.push_back(cc);
                    
                    double height_child = gtree->heights_gene[chi] > bpp.heights[ss] ? gtree->heights_gene[chi]:bpp.heights[ss];
                    double lentoCal=gtree->heights_gene[gs] - height_child;
                    if(lentoCal<0){lentoCal=1e-9;}
                    logTMs[cc] = bpp.getlogTMc(lentoCal, rates[rr], eigenvec, eigenval, eigeninv);
                }
               
                for(int g=0;g < len;g++)
                {
                    for(vector<int>::iterator cc=childs.begin(); cc<childs.end(); cc++)
                    {
                        int chi = gtree->children_gene[gs][*cc];
                        if(gtree->Tg[ss].at(chi)[g] >= bpp.num_base ) { //gtree->missing_gene[chi] || for speed
                            continue;
                        }                        
                        templambda[gs][g] += BPP::log_multi(logTMs[*cc], templambda[chi][g]);
                    }
                    //templambda[gs][g].col(2) = templambda[gs][g].col(1) + templambda[gs][g].col(0);
                }
            }
            
            if(sp < N)
            {
                for(map<int, int>::iterator it = gtree->parent_gene2[ss].begin(); it!= gtree->parent_gene2[ss].end(); it++)
                {
                    if(it->first != it->second) continue;
                    int chi = it->first;
                    if(gtree->missing_gene[chi]) continue;  //for speed
                    
                    double height_child = gtree->heights_gene[chi] > bpp.heights[ss] ? gtree->heights_gene[chi]:bpp.heights[ss];
                    double lentoCal=bpp.heights[sp] - height_child;
                    if(lentoCal<0){lentoCal=1e-9;}
                    mat log_TM = bpp.getlogTMc(lentoCal, rates[rr], eigenvec, eigenval, eigeninv);

                    for(int g = 0; g < len; g++)
                    {
                        if(gtree->Tg[ss].at(chi)[g] >= bpp.num_base) continue; //gtree->missing_gene[chi] || for speed
                        if(ss >= S || gtree->Tg[ss][chi][g] == -1)
                        {
                           // templambda[chi][g] = BPP::log_multi(log_TM, templambda[chi][g]);
                            emissions[rr] += BPP::log_exp_sum(log_TM.col(gtree->Tg[sp][chi][g]) + templambda[chi][g]);
                            //templambda[chi][g][gtree->Tg[sp].at(chi)[g]];
                        }else{
                            emissions[rr] += log_TM.at(gtree->Tg[ss][chi][g],gtree->Tg[sp][chi][g]);
                        }
                    }
                }
            }else{
                for(int g = 0; g < len; g++)
                {
                    if(gtree->missing_gene[ss] || gtree->Tg[ss].at(gtree->root)[g] >= bpp.num_base) continue;
                    emissions[rr] += templambda[gtree->root][g][gtree->Tg[ss].at(gtree->root)[g]];
                }
            }
            log_emission[ss][rr] =  emissions[rr];
        }   
    }   
}


void BPP_C::getUpdateNode(bool neut,vector<bool> & visited_init, BPP & bpp) //changedZ from small to large (bottom to top)
{

    for(vector<int>::iterator it = nodes.begin(); it < nodes.end();it++){
        int s = *it;
        if(missing[s])  continue;
        if(!neut && (Z[s]==0 ||  Z[s]==2)) continue;  // only get branches with Z == 1
        if(neut && (Z[s]==1 ||  Z[s]==0)) continue;   // only get branches with Z == 2
        int j = s;
        //j = bpp.parent[j];
        while(j!=N)
        {
            //j = bpp.parent[j];
            if(!visited_init[j]) break;
            visited_init[j] = 0;

            j = bpp.parent[j];
            assert(j>0);
        }
    }
}


double BPP_C::sample_rate(int indictor, int iter, int max_iter, vector<int> lens, int resZ, double old_rate, bool neut, vector<bool> visited, vec & loglik_old, BPP& bpp, int M, bool adaptive, double adaptive_factor)
{  //element number ,n or c, #MH steps
    //MH to sample rate

    double scale_adj;
    double r,cur_r = old_rate, proposal, MH_ratio,prior_a,prior_b;
    vec loglik_new(2);
    vec tmp_diag = zeros<mat>(bpp.num_base);

    if(neut) {prior_a = bpp.nprior_a; prior_b = bpp.nprior_b;} else {prior_a = bpp.cprior_a; prior_b = bpp.cprior_b;}

    if(visited[N-1]) {
   //     return(old_rate);
        // sample from prior
            cur_r = gsl_ran_gamma(RNG,prior_a,prior_b);
            if(neut && resZ!=0) {
                if(bpp.ropt == 1 && cur_r <  bpp.nlb ) //trace_c_rate[m]) 0.6
                {
                    cur_r = old_rate;
                }else if(bpp.ropt == 2 && cur_r < trace_c_rate[m])
                {
                    cur_r = old_rate;
                }
            }else if(!neut)
            {
                if(bpp.ropt == 1 && cur_r > bpp.cub) //trace_n_rate[m+1]) 1
                {
                    cur_r = old_rate;
                }else if(bpp.ropt == 2 && cur_r > trace_n_rate[m+1])
                {
                    cur_r = old_rate;
                }
            }
            else cur_r = old_rate;
    }else{
        int len = lens[indictor];
        vector<map<int, vector<mat>>> lambda_tmp = vector<map<int,vector<mat>>> (N, map<int,vector<mat>>()); // store orginial lambda

        for(int mm =0; mm<M;mm++)
        {
            r = gsl_rng_uniform(RNG);
            loglik_new = zeros<vec>(2);
            
            for(vector<int>::iterator it = nodes.begin(); it < nodes.end();it++)  //only copy visited[s]
            {
                int s = *it;
                if(s < S || visited[s])  continue;
                
                lambda_tmp[s] = gtree->lambda[s];
            }
            
            if(neut) {
                //double u = r*(prop_n - 1/prop_n ) + 1/prop_n; // generate uniform from (1/prop_n, prop_n);
                //proposal = cur_r *  u;
                proposal =gsl_ran_gamma(RNG,cur_r/prop_n,prop_n);
                if(bpp.ropt == 1 && proposal <  bpp.nlb ) //trace_c_rate[m]) 0.6
                {
                    continue;
                }else if(bpp.ropt == 2 && proposal < trace_c_rate[m])
                {
                     continue;
                }
            }else{
                double u = gsl_rng_uniform(RNG) *(prop_c - 1/prop_c ) + 1/prop_c; // generate uniform from (1/prop_n, prop_n);
                proposal = cur_r *  u;
                //proposal =gsl_ran_gamma(RNG,cur_r/prop_c,prop_c);
                if(bpp.ropt == 1 && proposal > bpp.cub) //trace_n_rate[m+1]) 1
                {
                    continue;
                }else if(bpp.ropt == 2 && proposal > trace_n_rate[m+1])
                {
                    continue;
                }

            }

          #pragma omp critical
          {

            if(proposal < 1e-12 || proposal > 1e8)
            {
                if(neut)
                {
                    cerr << "WARNING: sampling accelerated rate "<< proposal <<"  " << prop_n << " " <<cur_r <<" for element "<< CC <<"at iter " << m <<"  out of range!" <<endl;
                }else{
                    cerr << "WARNING: sampling conserved rate "<< proposal << "  " << prop_c <<"  "<< cur_r<<" for element "<< CC <<"at iter " << m << "  out of range!" <<endl;
                }
            }
          }

          if(proposal < 1e-12 || proposal > 1e8) {
          //    verbose=1;
              continue;
          }
           
           if(!neut)
           {
               gtree->Update_Tg(lens[indictor], visited, bpp, false, nodes, Z, trace_n_rate[m+1], proposal, log_pi, eigenvec, eigenval, eigeninv); //log_f_Xz(lens[1], visited, lambda_tmp, neut, proposal,bpp);  //P(X|Z), lambda not changed, lambda_tmp changed
           }else{
               gtree->Update_Tg(lens[indictor], visited, bpp, false, nodes, Z, proposal, trace_c_rate[m], log_pi, eigenvec, eigenval, eigeninv);
           }
            
            for(int g=0;g<lens[indictor];g++)
            {
                double temp = BPP::log_exp_sum(gtree->lambda[N-1][gtree->root][g].col(2) + log_pi);
                loglik_new[0] += temp;
                if(g < lens[1]) loglik_new[1] += temp;
            }
            
            if(neut) MH_ratio = loglik_new[indictor] -loglik_old[indictor] + log(gsl_ran_gamma_pdf(proposal,prior_a,prior_b)) - log(gsl_ran_gamma_pdf(cur_r,prior_a,prior_b)) + log(gsl_ran_gamma_pdf(cur_r,proposal/prop_n,prop_n)) - log(gsl_ran_gamma_pdf(proposal,cur_r/prop_n,prop_n)); //+ log(cur_r) - log(proposal);
            else  MH_ratio = loglik_new[indictor] -loglik_old[indictor] + log(gsl_ran_gamma_pdf(proposal,prior_a,prior_b)) - log(gsl_ran_gamma_pdf(cur_r,prior_a,prior_b)) + log(cur_r) - log(proposal);
            
            //cout << "proposal, cur:" << proposal <<", "<< cur_r << " log_lik_new: " << loglik_new[0] << ", " << loglik_old[0] << endl;

            if(log(r) > MH_ratio) // reject, restore
            {
                
                for(vector<int>::iterator it = nodes.begin(); it < nodes.end();it++)  //only copy visited[s]
                {
                    int s = *it;
                    if(s < S || visited[s])  continue;
                    
                    gtree->lambda[s] = lambda_tmp[s];
                }

            }else{  // update
                cur_r = proposal;
                for(int ii =0; ii <2; ii++) loglik_old[ii] = loglik_new[ii] ;
            }
        }
    }

    int new_burn=num_burn;
    //if(iter<max_iter) new_burn=floor(num_burn/2.0);
    if(m >=new_burn  &&  old_rate!= cur_r)
    {
        if(neut)
        {
            accept_n_rate +=1;
        }else{
            accept_c_rate +=1;

        }
    }

    //adaptive MCMC
    if (adaptive && m >=new_burn) {
        if ((m+1 - new_burn) % (adaptive_freq * num_thin) == 0) { //
            if(neut)
            {
                if((double) accept_n_rate/adaptive_freq > 0.44)
                    scale_adj = exp(min(adaptive_factor, 1 / sqrt((m+1) / (adaptive_freq * num_thin))));
                else if ((double) accept_n_rate/adaptive_freq < 0.23)
                    scale_adj = exp(-min(adaptive_factor, 1 / sqrt((m+1) / (adaptive_freq * num_thin))));
                else
                    scale_adj = 1;
                prop_n = prop_n * scale_adj;
                if(prop_n > 1) prop_n = 1;
                accept_n_rate = 0;
            }else{
                if((double) accept_c_rate/adaptive_freq > 0.44)
                    scale_adj = exp(min(adaptive_factor, 1 / sqrt((m+1) / (adaptive_freq * num_thin))));
                else if ((double) accept_c_rate/adaptive_freq < 0.23)
                    scale_adj = exp(-min(adaptive_factor, 1 / sqrt((m+1) / (adaptive_freq * num_thin))));
                else
                    scale_adj = 1;
                prop_c = pow(prop_c, scale_adj);
                if(prop_c < 1.00001) prop_c = 1.00001;
                
                accept_c_rate = 0;
            }
        }
    }

    return(cur_r);
}

vector<int> BPP_C::Update_Z(int len, BPP & bpp)  //whether cal prob_back again
{
    if(len >0 ) getEmission(len, bpp);
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end();it++) 
    {
        int s = *it;
        if(len >0 ) log_prob_back[s] = log_emission[s];
        else log_prob_back[s].fill(0);

        if(fixZ[s] ==1)
        {
            log_prob_back[s][2] = -INFINITY;
        }else if(fixZ[s] ==2)  // no use!
        {
            log_prob_back[s][0] = -INFINITY;
        }
    }

    //message passing from bottom to top
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s = *it;
        int p = bpp.parent[s];

        log_prob_back[p] += BPP::log_multi2(log_TM_Int[s], log_prob_back[s]);  //matrix * vector
    }

    vector<int> changedZ;

    //update Z from top to bottom Z[N] = 0
    int old;
    vec log_trans_p, trans_p;

    // sample root
    log_prob_back[N-1] += prior_z; //prior_z is only for root.

    vec prob = BPP::log_sample(log_prob_back[N-1]);

    unsigned int n[3];
    gsl_ran_multinomial(RNG, 3, 1,prob.memptr(),n);
    int nn;
    for(nn=0; nn<3;nn++)
    {
        if(n[nn]>0) break;
    }

    old = Z[N-1];

    Z[N-1] =nn;
    if(nn>1) {
        cout<< "Sample root Z error" <<endl;
    }

    if(old != nn)
    {
        changedZ.push_back(N-1);
    }

    for(vector<int>::iterator it = nodes.end()-2; it >= nodes.begin();it--)  // exclude root, doesn't include missing nodes
    {
        int s = *it;
        int p = bpp.parent[s];

        old = Z[s];
        //if(Z[p]!=2)
        //{
            log_trans_p =log_TM_Int[s].unsafe_col(Z[p]) + log_prob_back[s];
            trans_p = BPP::log_sample(log_trans_p);

            unsigned int n[3];
            gsl_ran_multinomial(RNG, 3, 1,trans_p.memptr(),n);
            int nn;
            for(nn=0; nn<3;nn++)
            {
                if(n[nn]>0) break;
            }
            Z[s] =nn;
            if(nn>2) {
                cout<< "Sample Z error" <<endl;
            }
//        }
//        else{
//            Z[s] =2;
//        }
        if(old != Z[s]) changedZ.push_back(s);
    }
    return(changedZ);
}

int BPP_C::sample_pi(int m, int indictor, vector<int> lens, int resZ, double piA_old, vec &loglik_old, double pi_delta , BPP &bpp, vector<double> obsCount, int resample_burnin){
    int pi_acc=0;
    vector<bool> visited_no = vector<bool>(N, false);
    vector<map<int, vector<mat> > > tmp_lambda = gtree->lambda;
    mat tmp_submat=submat;
    vec tmp_eigenval=eigenval;
    mat tmp_eigeninv=eigeninv;
    mat tmp_eigenvec=eigenvec;

    double tmpA;
    double log_p_prop;
    double log_p_cur;
    if(resample_burnin==1){
        tmpA = (double)gsl_ran_beta(RNG, prior_dir_param[0] + obsCount[0] + obsCount[3], prior_dir_param[1] + obsCount[1] + obsCount[2]) / 2;
    }else{
        tmpA = piA_old + (2.0 * gsl_rng_uniform(RNG) - 1) * pi_delta;  //is it correct? 0517, wrong when tmpA == piA_old == 0.05
        log_p_prop = -log(2.0 * pi_delta);
        if (tmpA < 0.05){
            //log_p_prop += log(0.05 - tmpA);
            log_p_prop +=log(0.05-(piA_old-pi_delta));
            tmpA = 0.05;
        }else if (tmpA > 0.45){
            //log_p_prop += log(tmpA - 0.45);
            log_p_prop += log((piA_old+pi_delta) - 0.45);
            tmpA = 0.45;
        }
        log_p_cur=-log(2.0*pi_delta);
        if(piA_old==0.05 && ((tmpA-pi_delta)<0.05)){
            log_p_cur+=log(0.05-(tmpA-pi_delta));
        }else if(piA_old==0.45 && ((tmpA+pi_delta)>0.45)){
            log_p_cur+=log(tmpA+pi_delta-0.45);
        }
    }
    
    if(tmpA == piA_old) return(1);
   
    pi[1] = pi[2] = 0.5 - tmpA;
    pi[0] = pi[3] =tmpA;
    log_pi=log(pi);

    getQmat(pi, inst_rate); //if not accepted, need to update back; and also below 

    vec loglik_new(2);
    loglik_new[0]=loglik_new[1]=0.0;

    gtree->Update_Tg(lens[indictor], visited_no, bpp, false, nodes, Z, trace_n_rate[m], trace_c_rate[m], log_pi, eigenvec, eigenval, eigeninv);
    for (int g = 0; g < lens[0]; g++){
        double temp = BPP::log_exp_sum(gtree->lambda[N - 1][gtree->root][g].col(2) + log_pi);
        if (g < lens[1])
            loglik_new[1] += temp;
        loglik_new[0] += temp;
    }
    
    if(resample_burnin==1){
        for (int ii = 0; ii < 2; ii++)
            loglik_old[ii] = loglik_new[ii];
        for (int i = 0; i < 4; i++){
            trace_pi[m + 1][i] = pi[i];
        }
        pi_acc=1;
        return(pi_acc);
    }

    double MH_ratio=loglik_new[indictor]-loglik_old[indictor];
    double temp1=log(gsl_ran_beta_pdf(2 * tmpA, prior_dir_param[0], prior_dir_param[1]))-log(gsl_ran_beta_pdf(2 * piA_old, prior_dir_param[0], prior_dir_param[1]));
    MH_ratio+=temp1;
    MH_ratio+=log_p_cur-log_p_prop;
    double r=gsl_rng_uniform(RNG);
    
    if(log(r)<MH_ratio){
        for (int ii = 0; ii < 2; ii++)
            loglik_old[ii] = loglik_new[ii];
        for (int i = 0; i < 4; i++){
            trace_pi[m + 1][i] = pi[i];
        }
        pi_acc=1;
    }else{
        for (int i = 0; i < 4; i++){
            pi[i]=trace_pi[m][i] ;
        }
        log_pi=log(pi);
        trace_pi[m+1]=trace_pi[m];
        submat=tmp_submat; 
        eigenval=tmp_eigenval;
        eigenvec=tmp_eigenvec;
        eigeninv=tmp_eigeninv;
        gtree->lambda=tmp_lambda;
    }
    return(pi_acc);
}


void BPP_C::printSptree(BPP& bpp, int l){
    cout << "print gt:" << endl;
    for (int s = 0; s < N; s++){
        cout << "gt=" << s << ", g_ht=" << gtree->heights_gene[s] << "\tpa=" << gtree->parent_gene[s] << "\tkid1=" << gtree->children_gene[s][0] << "\tkid2=" << gtree->children_gene[s][1] << endl;
    }
    cout << "print pg2 info:" << endl;
    for (int s = 0; s < N; s++){
        cout << "sp=" << s << ":";
        for (map<int, int>::iterator it = gtree->parent_gene2[s].begin(); it != gtree->parent_gene2[s].end(); it++){
            cout << "pg2: cl=" << it->first << " pa=" << it->second << endl;
        }
    }
    cout<<"print Tg info:"<<endl;
    for(int s=0; s<N; s++){
        cout <<"sp="<< s << ":";
        for(map<int,vector<int>>::iterator it = gtree->Tg[s].begin(); it != gtree->Tg[s].end(); it++){
            int temp_e=it->first;
            //cout<<"Tg element: "<<temp_e<<", size="<<gtree->Tg[s][temp_e].size()<<endl;
            cout<<"Tg element:"<<temp_e<<": ";
            for(int i=0; i<l; i++){
                if(i%5==0) cout<<"["<<i<<"]";
                cout<<gtree->Tg[s][temp_e][i];
            }
            cout<<endl;
        }
    }
}
