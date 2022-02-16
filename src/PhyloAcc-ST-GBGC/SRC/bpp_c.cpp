//
//  bpp_c.cpp
//  PhyloAcc_init2
//
//  Created by hzr on 4/19/16.
//  Copyright © 2016 hzr. All rights reserved.
//

#include "bpp_c.hpp"
#include "bpp.hpp"


mat BPP_C::getRate(mat Q0, double s, double b) // Qij: i -> j; ACGT
{
    if(s ==0 && b ==0) return(Q0);
    mat Q = Q0;
    mat B = zeros<mat>(Q.n_rows, Q.n_cols);
    B.at(0,1) = b; B.at(0,2) = b; B.at(3,1) = b; B.at(3,2) = b;
    B.at(1,0) = -b; B.at(1,3) = -b; B.at(2,0) = -b; B.at(2,3) = -b;
    B += s;
    B.transform( [](double val) { if(val == 0) return 1.0;
        return val/(1 - exp(-val));
    } );
    Q %= B;
    //Q %= (s + b*B + 0.0001)/(1 - exp(-(s+b*B + 0.0001)));
    colvec c = sum(Q,1);
    Q.diag() -=c;
    return(Q);
}



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



void BPP_C::getSubtree_missing(int ss, set<int> & upper, int child)  // traverse from root to children, include root
{
    
    if(upper.find(ss)==upper.end() || ss<S) return;
    
    int s = ss;
    while(upper.find(s)!=upper.end() && s>=S)
    {
        int* p = children2[s];
        //assert(p[0]!=-1 || p[1]!=-1);
        if(missing[p[0]] && !missing[p[1]])
        {
            parent2[p[0]] = -1;
            //root = s;
            distances2[p[1]] +=distances2[s];
            if(s == root)
            {
                root = p[1];
            }else{
                children2[parent2[ss]][child] = p[1];
                parent2[p[1]] = parent2[ss];
            }
            s = p[1];
        }
        else if(!missing[p[0]] && missing[p[1]])
        {
            parent2[p[1]] = -1;
            //root = s;
            distances2[p[0]] +=distances2[s];
            if(s == root)
            {
                root = p[0];
            }else{
                children2[parent2[ss]][child] = p[0];
                parent2[p[0]] = parent2[ss];
            }
            s = p[0];
        }
        else{
            
            getSubtree_missing(p[0], upper, 0);
            getSubtree_missing(p[1], upper, 1);
            return;
        }

    }
    
    

}


//void BPP_C::getEmission_ambig() //whether calculate null
//{
//    //vec tmp(4); vec tmp2(4); vec tmp3(4);
//    for(int s =0; s<S;s++) //N-1
//    {
//        if(missing[s]) continue;
//        for(int g=0;g <GG;g++)
//        {
//            if(Tg[g][s] !=-1) continue;
//            ambiguousS_null[g][s] = BPP::log_multi(log_cache_TM_null[s],lambda[g][s]);
//        }
//    }
//}


void BPP_C::initMCMC(int iter, BPP&bpp, int resZ)  // assign small prob from Z = 1 to missing
{
    
    isGB = vector<bool>(N,false);
    if(iter == 0)
    {
        //initalize Z, all Z are 1 except the root
        Z = vector<int>(N,-1);
        for(vector<int>::iterator it = nodes.begin();it<nodes.end();it++) Z[*it] = 0; //1
        Z[root] = 0;
    }else{
        
        Z = bpp.cur_Z[resZ][CC];
        //B = bpp.cur_B[resZ][CC];
        
    }
    
//    if(resZ >0) {
//        Z = bpp.Max_Z[0][CC];
//        isGB = bpp.Max_B[0][CC];
//    }
    
    //Z = bpp.uniqZ2[0];  //initilize randomly, follow Markov Chain
    
    fixZ = vector<int>(N,-1);
    for(vector<int>:: iterator it = upper_c.begin(); it < upper_c.end(); it++)
    {
        fixZ[*it] = 1;
    }
    
    fixZ[root] =0;  //0: 0,1: 0/1, -1: all
    //with prior, null model, all branches except outgroup are conserved, only have to sample Z  (0/1) of outgroup
    if(resZ==2)
    {
        //initiate log_prob_back not used
        
        for(int s=0; s<S;s++) //0/1
        {
            fixZ[s] = 1;
        }
        
        for(vector<int>::iterator it = bpp.target_species.begin(); it < bpp.target_species.end(); it++) //1/2
        {
            int s= *it;
            //restore restriction on 2
            fixZ[s] = -1;
            
        }
    }
    else if(resZ==0)
    {
        for(vector<int>::iterator it = nodes.begin();it<nodes.end()-1;it++)  //don't include root, subtree doesn't have the root edge
        {
            int s = *it;
            fixZ[s] = 1;
        }
//        for(int s=0; s<S;s++) //0/1
//        {
//            fixZ[s] = 1;
//        }
    }
    
    MaxLoglik = -INFINITY;
       

    trace_n_rate[0] = bpp.cur_nrate[CC]; //1;
    trace_c_rate[0] = bpp.cur_crate[CC]; //ratio; log_cache_TM_cons??
    trace_n_GB[0] = bpp.cur_nGB[CC]; //1;
    //trace_c_GB[0] = bpp.cur_cGB[CC]; //ratio; log_cache_TM_cons??
    
    // get Q matrix, need to modify for branch-wise B!!
    vector<mat> Qs(6);
    Qs[0] = bpp.submat; //bpp.cur_nGB[CC]);
    Qs[1] = getRate(bpp.submat, trace_c_rate[0], 0); //bpp.cur_nGB[CC]);
    Qs[2] = getRate(bpp.submat, trace_n_rate[0], 0); //bpp.cur_nGB[CC]);
    Qs[3] = getRate(bpp.submat, 0, trace_n_GB[0]); //bpp.cur_nGB[CC]);
    Qs[4] = getRate(bpp.submat, trace_c_rate[0], trace_n_GB[0]); //bpp.cur_nGB[CC]);
    Qs[5] = getRate(bpp.submat, trace_n_rate[0], trace_n_GB[0]); //bpp.cur_nGB[CC]);
    
    init_cache_TM(Qs);

    
    log_prob_back = vector<vec> (N, zeros<vec>(6));
    
    for(int s=S; s<N; s++)
        for(int g=0; g<GG; g++)
            lambda[g][s].fill(0);
    
   // getEmission_ambig();
    
    log_emission = vector<vector<double> >(N-1, vector<double>(6,0.0)); //Z,B
    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it ++) // For all nodes !!! ....only terminal nodes, S
    {
        int s = *it;
        int p = parent2[s];
        if(missing[p]) continue;
        
        if(missing[s]){ // don't take into acoount B at missing branches
            log_emission[s][2] = log(bpp.nconsToMis) ;//+ log(bpp.nGBprior_c);
            log_emission[s][1] = log(bpp.consToMis) ;//-INFINITY;
            log_emission[s][0] = log(bpp.nconsToMis);
            
            log_emission[s][5] = log(bpp.nconsToMis);// + log(1-bpp.nGBprior_c);
            log_emission[s][4] = log(bpp.consToMis);//-INFINITY;
            log_emission[s][3] = log(bpp.nconsToMis);
        }  // no initialize log_emission for not missing branches but will calculate in get emission
        
    }
   
     // get the probability of constraints on Z
    // constraint log_prob_back for fixZ
    for(int s=0; s<N;s++)
    {
        if(fixZ[s] == 1)
        {
            log_prob_back[s][2] = -INFINITY;
            log_prob_back[s][5] = -INFINITY;
            log_emission[s][2] = -INFINITY;
            log_emission[s][5] = -INFINITY;
        }
    }
    
    // message passing from bottom to top
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s = *it;
        int p = parent2[s];
        
        log_prob_back[p].subvec(0,2) += BPP::log_multi(log_TM_Int[s], log_prob_back[s]) - log(2);  //matrix * vector
        log_prob_back[p].subvec(3,5) += BPP::log_multi(log_TM_Int[s], log_prob_back[s]) - log(2);
        //cout << log_prob_back[p] <<endl;
    }
    logpZC = log_prob_back[root][0]; // prob of P(C)
    cout <<"logpZC: " << logpZC << endl;
    
}

// ofstream & outZ, no use!
void BPP_C::Gibbs(int iter, BPP &bpp, ofstream & outZ, string output_path,string output_path2,int resZ, bool UpR, bool UpHyper, double lrate_prop, double grate_prop)
{
    
    vector<int> changedZ;
    vector<bool> visited = vector<bool>(N,false);
    int accept = 0;
    int accept2 = 0;
    //std::fill(visited.begin()+S,visited.end(),false);
    
    //initial lambda
    for(int g=0; g<GG; g++)
        Update_Tg(g, visited, bpp, false);  // get history for each base and get lambda;
    
    // inital logp_Z
    logp_Z = prior_Z(Z) - logpZC;
//    for(vector<int>::iterator it = nodes.begin(); it < nodes.end(); it++)
//    {
//        if(missing[*it]) logp_Z += log_emission[*it][Z[*it]];// for missing, only Z no B
//    }
    
    // inital log_emission and log_prob_back including missing branches
   // initial_prob_back(bpp);
    double loglik_old;
    double elbo_old;
    bool move_r = true;
    
    for(m=0; m<(num_burn+num_mcmc); m++)
    {
        MonitorChain(m, move_r, accept,accept2, bpp, logp_Z, resZ);  //cout in each function
        //cout << m << ": " << log_lik_old << " " << logp_Z << endl;
        
        if(m == num_burn+num_mcmc-1) break;

        // get expectation of H_m|Z_m,r_m
        if(move_r)
        {
            for(vector<int>:: iterator it = nodes.begin(); it < nodes.end(); it++)//set all posterior_Tg to zero
            {
                posterior_Tg_next[*it].fill(0);
            }
            for(int g=0; g<GG; g++)
                Update_Tg(g, visited, bpp, true);  // get posterior_Tg_next for next iteration;
        }
        
        if(m==0) {
            posterior_Tg = posterior_Tg_next;
        }
        
        //elbo_old = ELBO(); // posterior_Tg, log_cache_TM
        elbo_old = initial_prob_back(bpp); //update log_emission and log_prob_back
        loglik_old  = log_prob_back[root][0];
        //if(resZ==2) cout << loglik_old << endl;
        
        //sample r & b, and then Z and B
        vector<mat> log_cache_TM_cons0 = log_cache_TM_cons;
        vector<mat> log_cache_TM_cons_B0 = log_cache_TM_cons_B;
        vector<mat> log_cache_TM_neut0 = log_cache_TM_neut;
        vector<mat> log_cache_TM_neut_B0 = log_cache_TM_neut_B;
        vector<mat> log_cache_TM_null_B0 = log_cache_TM_null_B;
        
        //vector<vector<double> > log_emission0 = log_emission; //Z,B, save log_emission if rejected
        //vector<vec> log_prob_back0 = log_prob_back; // save log_prob_back if rejected
        //double loglik_old0 = loglik_old; // save ELBO integrate out Z, B if rejected
        
        bool move_r2 = sample_rates(bpp,loglik_old, resZ, 3);  // 3 MH moves for r & b; loglik_old updated; log_emission updated; log_cache_TM updated
        move_r = update_rates(bpp, elbo_old, 3, false); // elbo_old, lambda, message updated; set adaptive = false
        //cout << move_r2;
        accept2 += move_r2;
        
        if(move_r){
            accept++;
           
        }else{
            //log_emission = log_emission0;
            //log_prob_back = log_prob_back0;
            //loglik_old = loglik_old0;
            log_cache_TM_cons = log_cache_TM_cons0;
            log_cache_TM_cons_B = log_cache_TM_cons_B0;
            log_cache_TM_neut = log_cache_TM_neut0;
            log_cache_TM_neut_B = log_cache_TM_neut_B0;
            log_cache_TM_null_B = log_cache_TM_null_B0;
        }
        
        posterior_Tg = posterior_Tg_next;
        
    }
    
    bpp.log_liks_Z[resZ][CC] = MaxLoglik;
    bpp.Max_Z[resZ][CC] = Max_Z;
    bpp.Max_B[resZ][CC] = Max_B;
    
//    if(resZ==0)
//    {
//        bpp.cur_crate[CC] = Max_c;
//        bpp.cur_nrate[CC] = Max_n;  //doesn't matter, initalize Z with no 2, random sample r_n from prior
//        bpp.cur_nGB[CC] = Max_b;
//    }
    
    if(UpHyper) {
        bpp.log_liks_curZ[CC] = trace_loglik[m - 1]; bpp.cur_Z[resZ][CC] = Z;
        // compute Z transition matrix for proposal
        double loss, gain;
        log_f_Z(Z, log_TM_Int, gain, loss);  // only do this for full model, no constraint
        bpp.MH_ratio_gain[CC] = -gain;
        bpp.MH_ratio_loss[CC] = -loss;
        
        
        for(vector<int>::iterator it = nodes.begin(); it <nodes.end(); it++)
        {
            int s = *it;
            double z = 1 - grate_prop;
            log_TM_Int[s](0,0) = log(z);
            log_TM_Int[s](1,0) = log(1 - z);
            //double y = exp(-lrate_prop *distances2[s]);
            double y = 1 - lrate_prop;
            log_TM_Int[s](1,1) = log(y);
            log_TM_Int[s](2,1) = log(1-y);
            //log_TM_Int[s] = log(TM_Int[s]);
        }
        
        
        
        for(vector<int>:: iterator it = upper_c.begin(); it!=upper_c.end();it++)
        {
            log_TM_Int[*it](1,1) = 0;
            log_TM_Int[*it](2,1) = log(0);
            //log_TM_Int[*it] = log(TM_Int[*it]);
        }
        
        
        log_f_Z(Z, log_TM_Int, gain, loss);
        bpp.MH_ratio_gain[CC] += gain;
        bpp.MH_ratio_loss[CC] += loss;
    
        bpp.cur_crate[CC] = trace_c_rate[m-1];
        bpp.cur_nrate[CC] = trace_n_rate[m-1];  //doesn't matter, initalize Z with no 2, random sample r_n from prior
        bpp.cur_nGB[CC] = trace_n_GB[m-1];
    }
    
    
    
    
   // if(verbose) //BPP::printZ(N, Max_Z,bpp.children);
     #pragma omp critical
    {
        if(verbose) {
            cout <<CC <<": " << Max_m << ", " << MaxLoglik << ", "<< trace_c_rate[Max_m] << ", " << trace_n_rate[Max_m] << ", " << trace_n_GB[Max_m] << endl;
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
            cout << endl << "B: ";
            for(int s=0;s<N;s++ )
            {
                if(Max_B[s])
                {
                    cout << s << ", ";
                }
            }
            
            cout << endl << endl;
        }
    }
    
    
    if(verbose) Output_sampling(iter, output_path2, bpp, resZ);  //has to be before Output_init!!

}













// only update probability of nodes above changedZ
void BPP_C::getUpdateNode(vector<int> changedZ, vector<bool> & visited_init) //changedZ from small to large (bottom to top)
{
    for(int i =0;i<N;i++)
        visited_init[i] = true;
    
    if(changedZ.size()==0) return;
    for(vector<int>::iterator it = changedZ.end()-1 ; it >= changedZ.begin(); --it){
        int j = *it;
        j = parent2[j];
        
        while(j!=N)
        {
            if(!visited_init[j]) break;
            visited_init[j] = 0;
            for(int g=0; g<GG; g++) lambda[g][j].zeros();  // only fathers of changedZ!
            
            j = parent2[j];
            assert(j!=-1);
            
        }
    }
    
}

void BPP_C::Update_Tg(int g, vector<bool> visited, BPP& bpp, bool tosample)  // impute base pair of each internal node (except missing nodes)
{
    
    if(!tosample)
    {
        // 1. sending the lambda msg from leaves bottom up through the network
        for(vector<int>::iterator it = internal_nodes.begin(); it <internal_nodes.end(); it++)  // didn't use leaves
        {
            int s = *it;
            if(visited[s] || missing[s])  continue; //  || Tg[g][s] == bpp.num_base
            int* p = children2[s];
            
            
            for(int cc=0;cc<2;cc++)
            {
                int chi = p[cc];
                if(missing[chi] ) continue; // || Tg[g][chi] == bpp.num_base
                switch (2 * Z[chi] + isGB[chi]) {
                    case 0:
                        message[g][chi] =  BPP::log_multi(log_cache_TM_null[chi],lambda[g][chi]);
                        break;
                    case 1:
                        message[g][chi] =  BPP::log_multi(log_cache_TM_null_B[chi],lambda[g][chi]);
                        break;
                    case 2:
                        message[g][chi] =  BPP::log_multi(log_cache_TM_cons[chi],lambda[g][chi]);
                        break;
                    case 3:
                        message[g][chi] =  BPP::log_multi(log_cache_TM_cons_B[chi],lambda[g][chi]);
                        break;
                    case 4:
                       message[g][chi] =  BPP::log_multi(log_cache_TM_neut[chi],lambda[g][chi]);
                        break;
                    case 5:
                       message[g][chi] =  BPP::log_multi(log_cache_TM_neut_B[chi],lambda[g][chi]);
                        break;
                        
                    default:
                        break;
                }
                lambda[g][s] += message[g][chi];
                
//                if(Z[chi] ==2) lambda[g][s] +=  BPP::log_multi(log_cache_TM_neut[chi],lambda[g][chi]);
//                else if (Z[chi] ==0) lambda[g][s] +=  BPP::log_multi(log_cache_TM_null[chi],lambda[g][chi]);
//                else lambda[g][s] +=  BPP::log_multi(log_cache_TM_cons[chi],lambda[g][chi]);
            }
            
        }
        
        // 2. processing the distribution of root species
        
        if(!visited[root]) lambda[g][root] += log_pi;
    }
    else
    {
        vec prob = BPP::log_sample_norm(lambda[g][root]);
        posterior_Tg_next[root].col(0) += prob;
        
        vector<vec> log_trans_p = vector<vec>(N, zeros<vec>(bpp.num_base));
        log_trans_p[root] = log(prob);
        
        mat trans_p(bpp.num_base, bpp.num_base);
        for(vector<int>::iterator it =nodes.end()-2; it>=nodes.begin(); it--)
        {
            int s = *it;
            int p = parent2[s];
            if(missing[s] ) continue; //|| Tg[g][s] == bpp.num_base
            
            switch (2 * Z[s] + isGB[s]) {
                case 0:
                    trans_p = log_cache_TM_null[s];  // row: s, col: p
                    //trans_p.each_col()+= log_trans_p[p] - message[g][s]; //m+
                    break;
                case 1:
                    trans_p = log_cache_TM_null_B[s];
                    //trans_p.each_col()+= log_trans_p[p] - message[g][s]; //m+
                    break;
                case 2:
                    trans_p = log_cache_TM_cons[s];
                    //trans_p.each_col()+= log_trans_p[p] - message[g][s]; //m+
                    break;
                case 3:
                    trans_p = log_cache_TM_cons_B[s];
                    //trans_p.each_col()+= log_trans_p[p] - message[g][s]; //m+
                    break;
                case 4:
                    trans_p = log_cache_TM_neut[s];
                    //trans_p.each_col()+= log_trans_p[p] - message[g][s]; //m+
                    break;
                case 5:
                    trans_p = log_cache_TM_neut_B[s];
                    //trans_p.each_col()+= log_trans_p[p] - message[g][s]; //m+
                    break;
                default:
                   break;
            }
            trans_p.each_row()+= (log_trans_p[p] - message[g][s]).t(); //m+
            trans_p.each_col() += lambda[g][s];
            log_trans_p[s] = BPP::log_exp_rowsum(trans_p);
            
            //cout << exp(trans_p) <<endl;
            //cout <<accu(exp(trans_p)) << endl;
            
            posterior_Tg_next[s] += exp(trans_p);  // row: s, col: p
//            double x =accu(exp(trans_p));
//            if(x < 0.99 || x > 1.01)
//            {
//                cout << x <<" ";
//                assert(x==1);
//            }
//
            
        }
    }
    
    
//    ctnutils::DispVector(Tg[g]);
//    cout<<endl;
    
}

void BPP_C::MonitorChain(int m, bool move_r, int & accept, int & accept2, BPP &bpp, const double add_loglik,const int resZ)  //calculate P(X|Z, TM(r) ),
{
    if(m==0)
    {
        trace_loglik[m] = 0;
        // calculate marginal likelihood,
        for(int g=0;g<GG;g++)
        {
            trace_loglik[m] += BPP::log_exp_sum(lambda[g][root]); // for X
            
        }
    }
    
    if(!move_r){
        trace_loglik[m] = trace_loglik[m-1];
        trace_full_loglik[m] = trace_full_loglik[m-1];
    }else{
        // add prior of Z and prior of r
        //double gain, loss;
        //log_f_Z(Z, log_TM_Int, gain, loss);
        trace_full_loglik[m] = trace_loglik[m] + add_loglik ;
        trace_full_loglik[m] += log(gsl_ran_gamma_pdf(-trace_c_rate[m],bpp.cprior_a,bpp.cprior_b));
        if(resZ !=0)
        {
            trace_full_loglik[m] += log(gsl_ran_gamma_pdf(trace_n_rate[m] + bpp.nprior_c,bpp.nprior_a,bpp.nprior_b));
        }
        
        int gBcount = 0;
        int ungBcount = 0;
        vector<int> mis_vs_Z = vector<int>(4,0); // record missing and Z==1, !mis & Z!=1,!mis & Z==1,  mis & Z!=1,  mis & Z==1,
        for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it ++)
        {
            int p = parent2[*it];
            if(missing[p]) continue;
            mis_vs_Z[2 *(missing[*it]) + (Z[*it] == 1)] += 1;
            if(missing[*it]) continue;
            if(isGB[*it] > 0)
            {
                gBcount ++;
                //trace_full_loglik[m] += log(gsl_ran_gamma_pdf(trace_n_GB[m],bpp.nGBprior_a,bpp.nGBprior_b)) + log(1 - bpp.nGBprior_c) ;
            }else{
                ungBcount ++;
                //trace_full_loglik[m] += log(bpp.nGBprior_c);
            }
        }
        trace_full_loglik[m] += log(1 - bpp.nconsToMis) * mis_vs_Z[0] + log(1 - bpp.consToMis) * mis_vs_Z[1] + log(bpp.nconsToMis) * mis_vs_Z[2] +log(bpp.consToMis) * mis_vs_Z[3];
        trace_full_loglik[m] += log(1 - bpp.nGBprior_c) * (ungBcount) + log(gsl_ran_gamma_pdf(trace_n_GB[m],bpp.nGBprior_a,bpp.nGBprior_b)) + log(bpp.nGBprior_c) * gBcount;
        
    }
    
    for(int s=0; s<N;s++)
    {
        trace_Z[m][s] = Z[s];
        trace_GB[m][s] = isGB[s];
    }
    
    
    if(m>=num_burn && trace_full_loglik[m]>MaxLoglik) //full
    {
        MaxLoglik = trace_full_loglik[m];
        Max_Z = Z;
        Max_m = m;
        Max_B = isGB;
        Max_b = trace_n_GB[m];
        Max_n = trace_n_rate[m];
        Max_c = trace_c_rate[m];
    }
    
    
    if(m%100==0 && verbose){
        cout <<CC <<": " << trace_loglik[m] <<", " << trace_full_loglik[m]<< ", " << trace_c_rate[m] << ", " << trace_n_rate[m] <<", " <<  prop_c <<", " << prop_n;
        cout << ", "  << trace_n_GB[m] <<", "  << prop_nGB << endl;
        cout <<"logp_Z: "<< logp_Z <<"log_prob_back: " <<log_prob_back[root][0] << endl;
        //cout << ", " << trace_c_GB[m] << ", " << trace_n_GB[m] <<", " <<  prop_cGB <<", " << prop_nGB << endl;
        cout << "accept ratio sb & ZB: " << (double)accept/100 << " " << (double)accept2/100 << endl;
        accept = 0; accept2 = 0;
        cout << "gBGC nodes: ";
        for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1; it++)
        {
            if(isGB[*it] > 0)
            {
                cout << *it << ", ";
            }
        }
       
        cout << endl;
        //cout << "changed Z: ";
        //ctnutils::DispVector(changedZ);
        //cout << "number of Z changed: " << changedZ.size();
        for(int s=0;s<N;s++ )
        {
            if(Z[s]!=1)
            {
                cout << s << ": " <<Z[s] << ", ";
            }
        }
        
        cout << endl;
    }

    
}

void BPP_C::getUpdateNode(bool neut, bool bg, vector<bool> & visited_init) //changedZ from small to large (bottom to top)
{
    
    //std::fill(visited_init, visited_init + N, 1);
    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end();it++){
        int s = *it;
        if(missing[s])  continue;
        if(!bg)
        {
            if(!neut && (Z[s]==0 ||  Z[s]==2)) continue;  // only get branches with Z == 1
            if(neut && (Z[s]==1 ||  Z[s]==0)) continue;   // only get branches with Z == 2
        }else{
            if(!isGB[s]) continue;
        }
        int j = s;
        j = parent2[j];
        
        while(j!=N)
        {
            if(!visited_init[j]) break;
            visited_init[j] = 0;
           
            j = parent2[j];
            assert(j>0);
        }
    }
}



//double BPP_C::log_f_Xz(vector<bool> visited, vector<vector<vec>>& lambda_tmp,vector<vector<vec>>& message_tmp, bool neut, bool gBC, double propos, BPP& bpp)  //calculate marginal likelihood P(X|Z,I)
//{
//    // 1. sending the lambda msg from leaves bottom up through the network
//    vec tmp_diag = zeros<mat>(bpp.num_base);
//    mat transition(bpp.num_base,bpp.num_base);
//    for(vector<int>::iterator it= internal_nodes.begin(); it < internal_nodes.end(); it ++ ) //int s=S; s<=root; s++)
//    {
//        int s = *it;
//        if(visited[s]) continue; //visited s includes missing s
//
//        int* p = children2[s];
//
//        for(int cc=0;cc<2;cc++)
//        {
//            int chi = p[cc];
//            assert(chi != -1);
//
//            if(missing[chi]) continue;
//
//            //Q = getRate(bpp.submat, propos, B);
//
//            if(gBC)
//            {
//                switch (2 * Z[chi] + isGB[chi]) {
//                    case 0:
//                        transition = log_cache_TM_null[chi]; break;
//                    case 1:
//                        transition = getTransition(chi, 1, propos, bpp);
//                        break;
//                    case 2:
//                        transition = log_cache_TM_cons[chi];
//                        break;
//                    case 3:
//                        transition = getTransition(chi, trace_c_rate[m+1], propos, bpp);
//                        break;
//                    case 4:
//                        transition = log_cache_TM_neut[chi]; break;
//                    case 5:
//                        transition = getTransition(chi, trace_n_rate[m+1], propos, bpp);
//                        break;
//                }
//
//            }else{
//                switch (2 * Z[chi] + isGB[chi]) {
//                    case 0:
//                        transition = log_cache_TM_null[chi]; break;
//                    case 1:
//                        transition = log_cache_TM_null_B[chi]; break;
//                    case 2:
//                        if(!neut) transition = getTransition(chi, propos, 0, bpp);
//                        else transition = log_cache_TM_cons[chi];
//                        break;
//                    case 3:
//                        if(!neut) transition = getTransition(chi, propos, trace_n_GB[m], bpp);
//                        else transition = log_cache_TM_cons_B[chi];
//                        break;
//                    case 4:
//                        if(neut) transition = getTransition(chi, propos, 0, bpp);
//                        else transition = log_cache_TM_neut[chi];
//                        break;
//                    case 5:
//                        if(neut) transition = getTransition(chi, propos, trace_n_GB[m], bpp);
//                        else transition = log_cache_TM_neut_B[chi];
//                        break;
//                }
//
//
////                if(!neut){  // sample cons rate
////                    if(Z[chi] ==1) {
////                        transition = getTransition(chi, propos, trace_n_GB[m+1], bpp);
////                    }
////                    else if(Z[chi] == 2) transition = log_cache_TM_neut[chi];
////                    else transition = log_cache_TM_null[chi];
////                }
////                else{
////                    if( Z[chi] ==2) {
////                       transition = getTransition(chi, propos, trace_n_GB[m+1], bpp);
////                    }
////                    else if (Z[chi] ==1) transition = log_cache_TM_cons[chi];
////                    else transition = log_cache_TM_null[chi];
////                }
//            }
//
//
//            for(int g=0;g<GG;g++)
//            {
//                message_tmp[g][chi] = BPP::log_multi(transition,lambda_tmp[g][chi]);
//                lambda_tmp[g][s] +=  message_tmp[g][chi]; //if(Tg[g][chi]!=bpp.num_base)
//
//            }
//
//        }
//    }
//
//
//
//
//    // 2. processing the distribution of root species
//    double result = 0;
//    for(int g=0;g<GG;g++)
//    {
//        if(!visited[root]) lambda_tmp[g][root] += log_pi;
//        result += BPP::log_exp_sum(lambda_tmp[g][root]);
//    }
//    //if(c==0) cout << lambda_tmp[0][N-1] <<endl;
//    return(result);
//}

double BPP_C::ELBO()
{
    double loglik = 0;
    mat trans_p;
   
    for(vector<int>::iterator it= nodes.begin(); it < nodes.end() - 1; it ++ ) //int s=S; s<=root; s++)
    {
        int s = *it;
        if(missing[s]) continue;
        
        switch (2 * Z[s] + isGB[s]) {
            case 0:
                trans_p = log_cache_TM_null[s];  // row: s, col: p
                break;
            case 1:
                trans_p = log_cache_TM_null_B[s];
                break;
            case 2:
                trans_p = log_cache_TM_cons[s];
                break;
            case 3:
                trans_p = log_cache_TM_cons_B[s];
                break;
            case 4:
                trans_p = log_cache_TM_neut[s];
                break;
            case 5:
                trans_p = log_cache_TM_neut_B[s];
                break;
            default:
                break;
        }
        trans_p %=posterior_Tg[s] ;
        loglik += accu(trans_p);
    }
    
    return(loglik);
}

// integrate out Z and B to sample rate
//update all parameters together by one MH step
bool BPP_C::sample_rates(BPP& bpp, double & loglik_old,  int resZ, int M)
{
    //double proposal, prior_a,prior_b; //cur_r = old_rate,
    bool move_r = 0;
    vector<vector<double> > log_emission_tmp = log_emission; //Z,B
    vector<vector<mat >> log_cache_TM_tmp = vector<vector<mat >>(5, vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base))); //3
    
    double cur_gB = trace_n_GB[m];
    double cur_rc = trace_c_rate[m];
    double cur_rn = trace_n_rate[m];
    
    
    //double elbo_old = ELBO(); // condition on Z, B and r, b
    
    /*************** sample r and b ***************/
    for(int mm =0; mm<M;mm++)
    {
        /*** 1. sample GB ***/
        double u = (gsl_rng_uniform(RNG) - 0.5) * prop_nGB * 2; //gsl_rng_uniform(RNG)*(prop_nGB - 1/prop_nGB ) + 1/prop_nGB;
        double new_gB = cur_gB + u;
        if(new_gB < 0.01) new_gB = cur_gB;
        //prior_a = bpp.nGBprior_a; prior_b = bpp.nGBprior_b;
        
        /*** 2. sample neut ***/
        //prior_a = bpp.nprior_a; prior_b = bpp.nprior_b;
        double new_neut = cur_rn;
        if(resZ != 0)
        {
            u = gsl_rng_uniform(RNG)*(prop_n - 1/prop_n ) + 1/prop_n; // generate uniform from (1/prop_n, prop_n);
            new_neut = (cur_rn + bpp.nprior_c) *  u - bpp.nprior_c;
            
            //u = (gsl_rng_uniform(RNG) - 0.5) * prop_n * 2;
            //new_neut = cur_rn + u;
            //if(new_neut < -bpp.nprior_c + 0.001 ) new_neut = cur_rn;
        }

        
        /*** 3. sample cons ***/
        //prior_a = bpp.cprior_a; prior_b = bpp.cprior_b;
        u = gsl_rng_uniform(RNG)*(prop_c - 1/prop_c ) + 1/prop_c;
        double new_cons = cur_rc * u;
        
        
        // update Q matrix
        vector<mat> Qs = vector<mat> (5);
        Qs[0] = getRate(bpp.submat,new_cons,0);
        Qs[1] = getRate(bpp.submat,new_neut,0);
        Qs[2] = getRate(bpp.submat,0,new_gB);
        Qs[3] = getRate(bpp.submat,new_cons,new_gB);
        Qs[4] = getRate(bpp.submat,new_neut,new_gB);
        
        integrate_ZB(bpp, 1, 0, Qs,log_emission_tmp, log_cache_TM_tmp);  // compute new transition matrix and integrate out Z & B at current paramters
        
        double MH_ratio = log_prob_back[root][0] - loglik_old + log(gsl_ran_gamma_pdf(new_gB,bpp.nGBprior_a,bpp.nGBprior_b)) - log(gsl_ran_gamma_pdf(cur_gB ,bpp.nGBprior_a,bpp.nGBprior_b)); //- log(new_gB) + log(cur_gB); //+ log(gsl_ran_gamma_pdf(cur_r,proposal/prop_nGB,prop_nGB)) - log(gsl_ran_gamma_pdf(proposal ,cur_r/prop_nGB,prop_nGB));
        // cout << "sample gBC: " << MH_ratio;
        
        if(resZ != 0)
        {
            MH_ratio += log(gsl_ran_gamma_pdf((new_neut + bpp.nprior_c),bpp.nprior_a,bpp.nprior_b)) - log(gsl_ran_gamma_pdf((cur_rn + bpp.nprior_c),bpp.nprior_a,bpp.nprior_b)) ;//- log(new_neut + bpp.nprior_c) + log(cur_rn + bpp.nprior_c); //+ log(gsl_ran_gamma_pdf((cur_r+bpp.nprior_c),(proposal + bpp.nprior_c)/prop_n,prop_n)) - log(gsl_ran_gamma_pdf((proposal + bpp.nprior_c),(cur_r+bpp.nprior_c)/prop_n,prop_n));
        }
        
        MH_ratio += log(gsl_ran_gamma_pdf(-new_cons,bpp.cprior_a,bpp.cprior_b)) - log(gsl_ran_gamma_pdf(-cur_rc,bpp.cprior_a,bpp.cprior_b)) - log(-new_cons) + log(-cur_rc); //+ log(gsl_ran_gamma_pdf(-cur_r,-proposal/prop_c,prop_c)) - log(gsl_ran_gamma_pdf(-proposal,-cur_r/prop_c,prop_c));
        
        //cout << "sample parameters : " << log_prob_back[root][0] - loglik_old << " " << MH_ratio << " " << new_gB << ", " << new_neut << ", " << new_cons << endl;
        if(log(gsl_rng_uniform(RNG)) < MH_ratio)
        {
            cur_gB = new_gB;
            cur_rn = new_neut;
            cur_rc = new_cons;
            loglik_old = log_prob_back[root][0];
            getEmission_update(0, 1, log_emission_tmp, log_cache_TM_tmp, bpp, 1); // change log_emission to log_emission_tmp, and update log_cache_TM
            move_r = true;
        }else{
            getEmission_update(0, 1, log_emission_tmp, log_cache_TM_tmp, bpp, 0); // restore log_emission_tmp
        }
        
    }
    
    trace_c_rate[m+1] = cur_rc;  //temporally
    trace_n_rate[m+1] = cur_rn;
    trace_n_GB[m+1] = cur_gB;
    return(move_r);
    
}
//void BPP_C::sample_rates(BPP& bpp, double & loglik_old, int M)
//{
//    double proposal, prior_a,prior_b; //cur_r = old_rate,
//
//    vector<vector<double> > log_emission_tmp = log_emission; //Z,B
//    vector<vector<mat >> log_cache_TM_tmp = vector<vector<mat >>(3, vector<mat >(N, zeros<mat> (bpp.num_base,bpp.num_base)));
//
//    double cur_gB = trace_n_GB[m];
//    double cur_rc = trace_c_rate[m];
//    double cur_rn = trace_n_rate[m];
//
//
//    //double elbo_old = ELBO(); // condition on Z, B and r, b
//
//    /*************** sample r and b ***************/
//    for(int mm =0; mm<M;mm++)
//    {
//        /*** 1. sample GB ***/
//        double u = gsl_rng_uniform(RNG)*(prop_nGB - 1/prop_nGB ) + 1/prop_nGB;
//        proposal = cur_gB * u;
//        prior_a = bpp.nGBprior_a; prior_b = bpp.nGBprior_b;
//
//        // update Q matrix
//        vector<mat> Qs = vector<mat> (3);
//        Qs[0] = getRate(bpp.submat,1,proposal);
//        Qs[1] = getRate(bpp.submat,cur_rc,proposal);
//        Qs[2] = getRate(bpp.submat,cur_rn,proposal);
//
//        integrate_ZB(bpp, 1, 0, Qs,log_emission_tmp, log_cache_TM_tmp);  // compute new transition matrix and integrate out Z & B at current paramters
//
//        double MH_ratio = log_prob_back[root][0] - loglik_old + log(gsl_ran_gamma_pdf(proposal,prior_a,prior_b)) - log(gsl_ran_gamma_pdf(cur_gB ,prior_a,prior_b)) - log(proposal) + log(cur_gB); //+ log(gsl_ran_gamma_pdf(cur_r,proposal/prop_nGB,prop_nGB)) - log(gsl_ran_gamma_pdf(proposal ,cur_r/prop_nGB,prop_nGB));
//      // cout << "sample gBC: " << MH_ratio;
//
//        if(log(gsl_rng_uniform(RNG)) < MH_ratio)
//        {
//            cur_gB = proposal;
//            loglik_old = log_prob_back[root][0];
//            getEmission_update(0, 1, log_emission_tmp, log_cache_TM_tmp, bpp, 1); // change log_emission to log_emission_tmp, and update log_cache_TM
//        }else{
//            getEmission_update(0, 1, log_emission_tmp, log_cache_TM_tmp, bpp, 0); // restore log_emission_tmp
//        }
//
//        /*** 2. sample neut ***/
//        prior_a = bpp.nprior_a; prior_b = bpp.nprior_b;
//        u = gsl_rng_uniform(RNG)*(prop_n - 1/prop_n ) + 1/prop_n; // generate uniform from (1/prop_n, prop_n);
//        proposal = (cur_rn + bpp.nprior_c) *  u - bpp.nprior_c;
////        if(bpp.ropt == 1 && proposal <  bpp.nlb ) //trace_c_rate[m]) 0.6
////        {
////            return(cur_rn);
////        }else if(bpp.ropt == 2 && proposal < trace_c_rate[m])
////        {
////            return(cur_rn);
////        }
//
//        Qs = vector<mat> (2);
//        Qs[0] = getRate(bpp.submat,proposal,0);
//        Qs[1] = getRate(bpp.submat,proposal,cur_gB);
//
//        integrate_ZB(bpp, 0, 1, Qs,log_emission_tmp, log_cache_TM_tmp);
//
//        MH_ratio = log_prob_back[root][0] -loglik_old + log(gsl_ran_gamma_pdf((proposal + bpp.nprior_c),prior_a,prior_b)) - log(gsl_ran_gamma_pdf((cur_rn + bpp.nprior_c),prior_a,prior_b)) - log(proposal + bpp.nprior_c) + log(cur_rn + bpp.nprior_c); //+ log(gsl_ran_gamma_pdf((cur_r+bpp.nprior_c),(proposal + bpp.nprior_c)/prop_n,prop_n)) - log(gsl_ran_gamma_pdf((proposal + bpp.nprior_c),(cur_r+bpp.nprior_c)/prop_n,prop_n));
//
//       // cout << "; nrate: " << MH_ratio;
//
//
//        if(log(gsl_rng_uniform(RNG)) < MH_ratio)
//        {
//            //move_r = true;
//            cur_rn = proposal;
//            loglik_old = log_prob_back[root][0];
//            getEmission_update(1, 0, log_emission_tmp, log_cache_TM_tmp, bpp, 1); // change log_emission to log_emission_tmp
//        }else{
//            getEmission_update(1, 0, log_emission_tmp, log_cache_TM_tmp, bpp, 0); // restore log_emission_tmp, need log_emission here!!
//        }
//
//        /*** 3. sample cons ***/
//        prior_a = bpp.cprior_a; prior_b = bpp.cprior_b;
//        u = gsl_rng_uniform(RNG)*(prop_c - 1/prop_c ) + 1/prop_c;
//        proposal = cur_rc * u;
////        if(bpp.ropt == 1 && proposal > bpp.cub) //trace_n_rate[m+1]) 1
////        {
////            return(cur_r);
////        }else if(bpp.ropt == 2 && proposal > trace_n_rate[m+1])
////        {
////            return(cur_r);
////        }
//
//        //Qs = vector<mat> (2);
//        Qs[0] = getRate(bpp.submat,proposal,0);
//        Qs[1] = getRate(bpp.submat,proposal,cur_gB);
//
//        integrate_ZB(bpp, 0, 0, Qs, log_emission_tmp, log_cache_TM_tmp);
//
//        MH_ratio = log_prob_back[root][0] -loglik_old + log(gsl_ran_gamma_pdf(-proposal,prior_a,prior_b)) - log(gsl_ran_gamma_pdf(-cur_rc,prior_a,prior_b)) - log(-proposal) + log(-cur_rc); //+ log(gsl_ran_gamma_pdf(-cur_r,-proposal/prop_c,prop_c)) - log(gsl_ran_gamma_pdf(-proposal,-cur_r/prop_c,prop_c));
//
//       // cout << "; crate: " << MH_ratio << endl;
//
//        if(log(gsl_rng_uniform(RNG)) < MH_ratio)
//        {
//            //move_r = true;
//            cur_rc = proposal;
//            loglik_old = log_prob_back[root][0];
//            getEmission_update(0, 0, log_emission_tmp, log_cache_TM_tmp, bpp, 1); // change log_emission (will be used to restore) to log_emission_tmp
//        }else{
//            getEmission_update(0, 0, log_emission_tmp, log_cache_TM_tmp, bpp, 0); // restore log_emission_tmp
//        }
//    }
//
//    trace_c_rate[m+1] = cur_rc;  //temporally
//    trace_n_rate[m+1] = cur_rn;
//    trace_n_GB[m+1] = cur_gB;
//
//}

bool BPP_C::update_rates(BPP& bpp, double & elbo_old, int M, bool adaptive, double adaptive_factor)
{
    /*************** sample Z and B ***************/
    vector<int> changedZ = Update_ZB();  // no use log_cache_TM_tmp and Qs
    bool move_r = false;
    
    // compute true loglik integrate H given Z/B
    vector<bool> visited = vector<bool>(N,false);
    double loglik_new = 0;
    vector<vector<vec>> lambda0 = lambda;  // save for restore
    vector<vector<vec>> message0 = message; // save for restore
    
    for(int s=S; s<N; s++)
        for(int g=0; g<GG; g++)
            lambda[g][s].fill(0);
    
    for(int g=0; g<GG; g++)
    {
        Update_Tg(g, visited, bpp, false);  // get history for each base and get lambda;
        loglik_new += BPP::log_exp_sum(lambda[g][root]); // for X
    }
    double elbo_new = ELBO();
    double MH_ratio = loglik_new - trace_loglik[m] + elbo_old - elbo_new;
        
    //cout << "MH all: " << changedZ.size() << "; " <<  elbo_new  -  elbo_old << "; " << loglik_new  -  trace_loglik[m] << ";" << MH_ratio << endl;
        // cout <<"rates: " <<  cur_r << "; " << proposal <<"; " <<  loglik_delta << ";" <<loglik_new -loglik_old<< ";" << MH_ratio << "; " << time_approx*1000 << "; " << time_exact*1000 << "; " << endl;
    
    if(log(gsl_rng_uniform(RNG)) < MH_ratio)
    {
//        trace_c_rate[m+1] = cur_rc;
//        trace_n_rate[m+1] = cur_rn;
//        trace_n_GB[m+1] = cur_gB;
        trace_loglik[m+1] = loglik_new;
        if(changedZ.size() >0) logp_Z = prior_Z(Z) - logpZC;
        ////elbo_old = elbo_new;
        move_r = true;

    }else{
        //log_emission = log_emission0;
        //log_prob_back = log_prob_back0;
        //loglik_old = loglik_old0;
        trace_c_rate[m+1] = trace_c_rate[m];
        trace_n_rate[m+1] = trace_n_rate[m];
        trace_n_GB[m+1] = trace_n_GB[m];
       
        lambda = lambda0;
        message = message0;
        Z = trace_Z[m];
        isGB = trace_GB[m];
    }
    
    
    //adaptive MCMC
    if (adaptive) {
        adaptive_MCMC(adaptive_factor, 1); //M
    }
   
    return(move_r) ;
}

void BPP_C:: adaptive_MCMC(double adaptive_factor, int M)
{
    double scale_adj;
    accept_c_rate += (trace_c_rate[m+1] !=  trace_c_rate[m]);
    accept_n_rate += (trace_n_rate[m+1] !=  trace_n_rate[m]);
    accept_n_GB += (trace_n_GB[m+1]!=  trace_n_GB[m]);
    
    if ((m+1) % adaptive_freq == 0) {
        if(verbose)
        {
            cout << "accept ratios, c: " << (double)accept_c_rate/adaptive_freq <<" n: " << (double)accept_n_rate/adaptive_freq << " gb: " <<(double)accept_n_GB/adaptive_freq << endl << endl;
        }
        
        if((double) accept_n_rate/adaptive_freq > 1 - pow((1 -0.44), M))
            scale_adj = exp(min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
        else if ((double) accept_n_rate/adaptive_freq < 1 - pow((1 -0.23), M))
            scale_adj = exp(-min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
        else
            scale_adj = 1;
        prop_n = prop_n * scale_adj;
        if(prop_n < 1.01) prop_n = 1.01;
        accept_n_rate = 0;
        
        if((double) accept_c_rate/adaptive_freq > 1 - pow((1 -0.44), M))
            scale_adj = exp(min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
        else if ((double) accept_c_rate/adaptive_freq < 1 - pow((1 -0.23), M))
            scale_adj = exp(-min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
        else
            scale_adj = 1;
        prop_c = prop_c * scale_adj;
        if(prop_c < 1.01) prop_c = 1.01;
        accept_c_rate = 0;
        
        if((double) accept_n_GB/adaptive_freq > 1 - pow((1 -0.44), M))
            scale_adj = exp(min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
        else if ((double) accept_n_GB/adaptive_freq < 1 - pow((1 -0.23), M))
            scale_adj = exp(-min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
        else
            scale_adj = 1;
        prop_nGB = prop_nGB * scale_adj;
        if(prop_nGB < 1.01) prop_nGB = 1.01;
        accept_n_GB = 0;
    }
}


//double BPP_C::sample_rate(int resZ, double old_rate, bool neut, vector<bool> visited, double & loglik_old, BPP& bpp, int M, bool adaptive, double adaptive_factor)
//{  //element number ,n or c, #MH steps
//    //MH to sample rate
//    
//    double scale_adj;
//    double r,cur_r = old_rate, proposal, MH_ratio,loglik_new ,prior_a,prior_b;
//    
//    if(neut) {
//        prior_a = bpp.nprior_a; prior_b = bpp.nprior_b;
//        
//    } else {
//        prior_a = bpp.cprior_a; prior_b = bpp.cprior_b;
//    }
//    
//    if(visited[root]) {
//        // sample from prior
//        if(resZ!=0)
//        {
//            old_rate = gsl_ran_gamma(RNG,prior_a,prior_b);
//            if(neut) old_rate -= bpp.nprior_c;
//            else old_rate = -old_rate;
//        }
//        return(old_rate);
////        return(ratio);  // originally for acceleration
////        if(neut)
////        {
////            return(1);
////        }else{
////            return(ratio);
////        }
//    }
//    
//    
//    vec tmp_diag = zeros<mat>(bpp.num_base);
//    vector<vector<vec>> lambda_tmp = vector<vector<vec>> (GG, vector<vec>(N,zeros<vec>(bpp.num_base)));
//    vector<vector<vec>> message_tmp = vector<vector<vec>> (GG, vector<vec>(N,zeros<vec>(bpp.num_base)));
//    
//    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)  // only copy visited[s]
//    {
//        int s = *it;
//        if(!visited[s]) continue;
//        for(int g=0;g<GG;g++) lambda_tmp[g][s] = lambda[g][s];
//    }
//    
//    
//    
//    for(int mm =0; mm<M;mm++)
//    {
//        r = gsl_rng_uniform(RNG);
//        
//        if(neut) {
//            //proposal = gsl_ran_gamma(RNG,(cur_r + bpp.nprior_c)/prop_n,prop_n) - bpp.nprior_c;
//            double u = gsl_rng_uniform(RNG)*(prop_n - 1/prop_n ) + 1/prop_n; // generate uniform from (1/prop_n, prop_n);
//            proposal = (cur_r + bpp.nprior_c) *  u - bpp.nprior_c;
//            if(bpp.ropt == 1 && proposal <  bpp.nlb ) //trace_c_rate[m]) 0.6
//            {
//                return(cur_r);
//            }else if(bpp.ropt == 2 && proposal < trace_c_rate[m])
//            {
//                 return(cur_r);
//            }
//        }
//        else{
//            double u = gsl_rng_uniform(RNG)*(prop_c - 1/prop_c ) + 1/prop_c;
//            proposal = cur_r * u;
//            if(bpp.ropt == 1 && proposal > bpp.cub) //trace_n_rate[m+1]) 1
//            {
//                return(cur_r);
//            }else if(bpp.ropt == 2 && proposal > trace_n_rate[m+1])
//            {
//                return(cur_r);
//            }
//        }
//        time_t start = time(NULL);
//        loglik_new = log_f_Xz(visited, lambda_tmp, message_tmp, neut, 0, proposal,bpp);  //P(X|Z), lambda not changed, lambda_tmp changed
//        double time_exact = (time(NULL)-start)*1000;
//         // compare with loglik computed using ELBO
//        double loglik_delta = 0;
//         mat transition;
//        start = time(NULL);
//        for(vector<int>::iterator it= nodes.begin(); it < nodes.end() - 1; it ++ ) //int s=S; s<=root; s++)
//        {
//            int s = *it;
//           
//            //compute new transition matrix
//            if(neut && Z[s] == 2)
//            {
//                transition = getTransition(s, proposal, isGB[s]* trace_n_GB[m], bpp);
//                transition -= isGB[s]? log_cache_TM_neut_B[s]:  log_cache_TM_neut[s];
//            }else if(!neut && Z[s] == 1)
//            {
//                transition = getTransition(s, proposal, isGB[s]* trace_n_GB[m], bpp) ;
//                transition -= isGB[s]? log_cache_TM_cons_B[s]:  log_cache_TM_cons[s];;
//            }else{
//                continue;
//            }
//            transition %=posterior_Tg[s] ;
//            loglik_delta += accu(transition);
//        }
//         double time_approx = (time(NULL)-start);
//
//        if(neut) MH_ratio = loglik_new -loglik_old + log(gsl_ran_gamma_pdf((proposal + bpp.nprior_c),prior_a,prior_b)) - log(gsl_ran_gamma_pdf((cur_r + bpp.nprior_c),prior_a,prior_b)) - log(proposal + bpp.nprior_c) + log(cur_r + bpp.nGBprior_c); //+ log(gsl_ran_gamma_pdf((cur_r+bpp.nprior_c),(proposal + bpp.nprior_c)/prop_n,prop_n)) - log(gsl_ran_gamma_pdf((proposal + bpp.nprior_c),(cur_r+bpp.nprior_c)/prop_n,prop_n));
//        else  MH_ratio = loglik_new -loglik_old + log(gsl_ran_gamma_pdf(-proposal,prior_a,prior_b)) - log(gsl_ran_gamma_pdf(-cur_r,prior_a,prior_b)) - log(-proposal) + log(-cur_r); //+ log(gsl_ran_gamma_pdf(-cur_r,-proposal/prop_c,prop_c)) - log(gsl_ran_gamma_pdf(-proposal,-cur_r/prop_c,prop_c));
//        
//      
//        //cout << cur_r << "; " << proposal <<"; " << loglik_old  << "; " << loglik_new << ";" <<loglik_new -loglik_old<< ";" << MH_ratio << endl;
//       // cout <<"rates: " <<  cur_r << "; " << proposal <<"; " <<  loglik_delta << ";" <<loglik_new -loglik_old<< ";" << MH_ratio << "; " << time_approx*1000 << "; " << time_exact*1000 << "; " << endl;
//       
//        
//        if(log(r) < MH_ratio)
//        {
//            cur_r = proposal;
//            loglik_old = loglik_new;
//            //getUpdate some lambdas
//            for(vector<int>::iterator it= internal_nodes.begin(); it < internal_nodes.end(); it ++ )  // int s=S; s<=root;s++)
//            {
//                int s = *it;
//                if(visited[s]) continue; // || parent2[s] == -1
//                for(int g=0;g<GG;g++)
//                {                    
//                    lambda[g][s] = lambda_tmp[g][s];
//                    
//                }
//                
//                int* p = children2[s];
//                for(int cc=0;cc<2;cc++)
//                {
//                    int chi = p[cc];
//                    assert(chi != -1);
//                    
//                    if(missing[chi]) continue;
//                    for(int g=0;g<GG;g++)
//                    {
//                        message[g][chi] = message_tmp[g][chi];
//                        
//                    }
//                }
//            }
//            
//            
//        }
//        
//        // reset lambda_tmp to zero
//        if(mm<M-1)
//        {
//            for(vector<int>::iterator it= internal_nodes.begin(); it < internal_nodes.end(); it ++ )//int s=S; s<=root;s++)
//            {
//                int s = *it;
//                if(visited[s]) continue; //|| parent2[s] == -1
//                for(int g=0;g<GG;g++)
//                {
//                    lambda_tmp[g][s].fill(0);
//                }
//            }
//        }
//        
//    }
//    
//    if(old_rate!= cur_r)
//    {
//        // update Q matrix
//        vector<mat> Qs = vector<mat> (2);
//        Qs[0] = getRate(bpp.submat,cur_r,0);
//        Qs[1] = getRate(bpp.submat,cur_r,trace_n_GB[m]);
//       
//        update_cache_TM(Qs, 0 ,neut);
//        
////        for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
////        {
////            int s = *it;
////            //if(s==bpp.moveroot) continue; //not update rate for outgroup! s==42||, conserved rate still changing
////            if(neut) log_cache_TM_neut[s] = getTransition(s, cur_r, trace_n_GB[m+1], bpp); //transpose Q
////            else log_cache_TM_cons[s] = getTransition(s, cur_r, trace_n_GB[m+1], bpp);
////        }
//        
//        //getEmission_ambig(false);
//        if(neut)
//        {
//            accept_n_rate +=1;
//        }else{
//            accept_c_rate +=1;
//            
//        }
//        
//    }
//    
//    //adaptive MCMC
//    if (adaptive) {
//        if ((m+1) % adaptive_freq == 0) {
//            if(neut)
//            {
//                if((double) accept_n_rate/adaptive_freq > 0.44)
//                    scale_adj = exp(min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
//                else if ((double) accept_n_rate/adaptive_freq < 0.23)
//                    scale_adj = exp(-min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
//                else
//                    scale_adj = 1;
//                prop_n = prop_n * scale_adj;
//                if(prop_n < 1.01) prop_n = 1.01;
//                accept_n_rate = 0;
//            }else{
//                if((double) accept_c_rate/adaptive_freq > 0.44)
//                    scale_adj = exp(min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
//                else if ((double) accept_c_rate/adaptive_freq < 0.23)
//                    scale_adj = exp(-min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
//                else
//                    scale_adj = 1;
//                prop_c = prop_c * scale_adj;
//                if(prop_c < 1.01) prop_c = 1.01;
//                accept_c_rate = 0;
//            }
//        }
//    }
//    
//
//    
//    return(cur_r);
//}

// sample B per branch
mat BPP_C::getTransition(int s, double sel, double bias, BPP & bpp)
{
    mat Q;
    Q = getRate(bpp.submat, sel, bias);

    cx_mat bvec;
    cx_vec aval;
    eig_gen(aval, bvec, Q);
    vec eigenval = conv_to<mat>::from(aval);
    mat eigenvec = conv_to<mat>::from(bvec).t();
    mat eigeninv;
    try{
        eigeninv = inv(eigenvec);
    }catch(std::runtime_error & e){
        cout << "get transition error: " << Q <<endl;
    }
    
    mat x = eigenvec;
    vec tmp_diag  = exp(eigenval*distances2[s]);
    x.each_col()%=tmp_diag;
    return(log(eigeninv * x));

}

//double BPP_C::update_B(int s, double old_rate, BPP& bpp, int M, bool adaptive, double adaptive_factor)
//{  //element number ,n or c, #MH steps
//    //MH to sample rate
//
//    double scale_adj;
//    double r,cur_r = old_rate, proposal, MH_ratio,loglik_new ,prior_a,prior_b;
//    prior_a = bpp.nGBprior_a; prior_b = bpp.nGBprior_b;
//
//    // inital estimate of Bs;
//    int p = parent2[s];
//    mat subs = zeros<mat>(bpp.num_base, bpp.num_base);  // A, C, G, T
//
//    //compute number of different substitions
//    for(int g = 0; g <GG; g++)
//    {
//        if(Tg[g][s] == bpp.num_base || Tg[g][s] == -1) continue;
//        subs.at(Tg[g][p], Tg[g][s]) += 1;
//    }
//
//    if(m == 0 ) // set initial value of B
//    {
//        vec bcomp = sum(subs, 1);   // number of ACGT at parent;
//        double ratio = (subs.at(0,1) + subs.at(0,2) + subs.at(3,1) + subs.at(3,2) + 1)/ (subs.at(1,0) + subs.at(2,0) + subs.at(1,3) + subs.at(2,3) + 1); // ratio of nAT->CG/nCG->AT
//        double ratio0 = bcomp.at(0) * (bpp.submat.at(0,1) + bpp.submat.at(0,2)) + bcomp.at(3) * (bpp.submat.at(3,1) + bpp.submat.at(3,2));
//        ratio0/= (bcomp.at(1) * (bpp.submat.at(1,0) + bpp.submat.at(1,3)) + bcomp.at(2) * (bpp.submat.at(2,0) + bpp.submat.at(2,3)));// ratio of nAT->CG/nCG->AT with no bias;
//        proposal = log(ratio) - log(ratio0);
//        cout << s << ": " << proposal << endl;
//        cout << subs << endl;
//        if(proposal <= 1e-4)
//        {
//            proposal = 0;
//            //isGB[s]  = 0;
//        }
//        cur_r = proposal;
//
//    }else
//    {
//        double loglik_old, sel;
//        switch (Z[s])
//        {
//            case 1:
//                sel = trace_c_rate[m];
//                loglik_old = accu(log_cache_TM_cons[s] % subs);
//                break;
//            case 2:
//                sel = trace_n_rate[m];
//                loglik_old = accu(log_cache_TM_neut[s] % subs);
//                break;
//            default:
//                loglik_old = accu(log_cache_TM_null[s] % subs);
//                sel = 1;
//        }
//
//
//        //mat transition_prop;
//        for(int mm =0; mm<M;mm++)
//        {
//            r = gsl_rng_uniform(RNG);
//            //proposal = gsl_ran_gamma(RNG,cur_r /prop_nGB[s],prop_nGB[s]);
//            //MH_ratio = 0;
//           if(cur_r > 0)
//           {
//               if(r > prop_nGB2[s]){
//                   proposal = gsl_ran_gaussian(RNG, prop_nGB[s]) + log(cur_r);  // log-normal
//                   MH_ratio = -log(cur_r) + proposal ; // proposal already taken log
//                   proposal = exp(proposal);
//                   MH_ratio +=  log(gsl_ran_gamma_pdf(proposal,prior_a,prior_b)) - log(gsl_ran_gamma_pdf(cur_r ,prior_a,prior_b));
//               }
//               else {
//                   proposal = 0;
//                   MH_ratio = -log(prop_nGB2[s]) + log(gsl_ran_gamma_pdf(cur_r ,prior_a,prior_b));
//                   MH_ratio -=  log(gsl_ran_gamma_pdf(cur_r,prior_a,prior_b)) + log(1 - bpp.nGBprior_c) - log(bpp.nGBprior_c);
//               }
//
//           }else{
//               proposal = gsl_ran_gamma(RNG,prior_a,prior_b); // sample from prior
//               MH_ratio = log(prop_nGB2[s]) - log(gsl_ran_gamma_pdf(proposal ,prior_a,prior_b));
//               MH_ratio +=  log(gsl_ran_gamma_pdf(proposal,prior_a,prior_b)) + log(1 - bpp.nGBprior_c) - log(bpp.nGBprior_c);
//           }
//
//            //if(proposal < 1e-4) continue;
//            mat transition = getTransition(s,sel, proposal, bpp);
//
//            loglik_new = accu(transition % subs);
//
//            MH_ratio += loglik_new -loglik_old;// + log(gsl_ran_gamma_pdf(proposal,prior_a,prior_b)) - log(gsl_ran_gamma_pdf(cur_r ,prior_a,prior_b));
//            //log(gsl_ran_gamma_pdf(cur_r,proposal/prop_nGB[s],prop_nGB[s])) - log(gsl_ran_gamma_pdf(proposal ,cur_r/prop_nGB[s],prop_nGB[s]));
//
//           //cout << cur_r << "; " << proposal <<"; " << loglik_old  << "; " << loglik_new << ";" <<loglik_new -loglik_old<< ";" << MH_ratio << endl;
//
//             r = gsl_rng_uniform(RNG);
//            if(log(r) < MH_ratio)
//            {
//                cur_r = proposal;
//                loglik_old = loglik_new;
//                //transition_prop = transition;
//
//            }
//
//        }
//    }
//
//    if(old_rate!= cur_r)
//    {
////        if(cur_r < 0.01)
////        {
////            isGB[s] = 0;
////            cur_r = 0;
////        }
//        // update Q matrix
//        log_cache_TM_cons[s] = getTransition(s,trace_c_rate[m], cur_r, bpp);
//        log_cache_TM_neut[s] = getTransition(s,trace_n_rate[m], cur_r, bpp);
//        log_cache_TM_null[s] = getTransition(s,1, cur_r, bpp);
//        accept_n_GB[s] +=1;
//
//    }
//
//    //adaptive MCMC
//    if (adaptive) {
//        if ((m+1) % adaptive_freq == 0) {
//                if((double) accept_n_GB[s]/adaptive_freq > 0.44)
//                    scale_adj = exp(min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
//                else if ((double) accept_n_GB[s]/adaptive_freq < 0.23)
//                    scale_adj = exp(-min(adaptive_factor, 1 / sqrt((m+1) / adaptive_freq)));
//                else
//                    scale_adj = 1;
//                prop_nGB[s] = prop_nGB[s] * scale_adj;
//                accept_n_GB[s] = 0;
//        }
//    }
//
//    return(cur_r);
//}

void BPP_C:: init_cache_TM(vector<mat> & Qs)
{
    size_t sz = Qs.size();
    vector<vec> eigenvals = vector<vec>(sz);
    vector<mat> eigenvecs =  vector<mat>(sz);
    vector<mat> eigeninvs =  vector<mat>(sz);
    
    vec tmp_diag;
    cx_mat bvec;
    cx_vec aval;
    for(int i =0; i< sz; i++)
    {
        eig_gen(aval, bvec, Qs[i]);
        eigenvals[i] = conv_to<mat>::from(aval);
        eigenvecs[i] = conv_to<mat>::from(bvec).t();
        try{
            eigeninvs[i] = inv(eigenvecs[i]);
        }catch(std::runtime_error & e){
            cout << "get transition error: " << Qs[i] <<endl;
        }
        
    }
    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s = *it;
        if(missing[s]) continue;
        //if(s==bpp.moveroot) continue; //not update rate for outgroup! s==42||, conserved rate still changing
        for(int i =0; i< sz; i++)
        {
            tmp_diag  = exp(eigenvals[i]*distances2[s]);
            mat x =eigenvecs[i];
            x.each_col()%=tmp_diag;
                switch (i) {
                    case 0:
                        log_cache_TM_null[s] = log(eigeninvs[i] * x);
                        break;
                    case 1:
                        log_cache_TM_cons[s] = log(eigeninvs[i] * x);
                        break;
                    case 2:
                        if(fixZ[s] != 1) log_cache_TM_neut[s] = log(eigeninvs[i] * x);
                        break;
                    case 3:
                        log_cache_TM_null_B[s] = log(eigeninvs[i] * x);
                        break;
                    case 4:
                        log_cache_TM_cons_B[s] = log(eigeninvs[i] * x);
                        break;
                    case 5:
                        if(fixZ[s] != 1) log_cache_TM_neut_B[s] = log(eigeninvs[i] * x);
                        break;
                    default:
                        break;
                }

        }
    }
}

void BPP_C:: update_cache_TM(vector<mat> & Qs, vector<vector<mat >>& log_cache_TM_tmp)
{
    size_t sz = Qs.size();
    vector<vec> eigenvals = vector<vec>(sz);
    vector<mat> eigenvecs =  vector<mat>(sz);
    vector<mat> eigeninvs =  vector<mat>(sz);
    
    vec tmp_diag;
    cx_mat bvec;
    cx_vec aval;
    for(int i =0; i< sz; i++)
    {
        eig_gen(aval, bvec, Qs[i]);
        eigenvals[i] = conv_to<mat>::from(aval);
        eigenvecs[i] = conv_to<mat>::from(bvec).t();
        try{
             eigeninvs[i] = inv(eigenvecs[i]);
        }catch(std::runtime_error & e){
            cout << "get transition error: " << Qs[i] <<endl;
        }
        
    }
    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s = *it;
        if(missing[s]) continue;

        for(int i =0; i< sz; i++)
        {
            //if((fixZ[s] == 1) & gB & i == sz -1 ) continue;
            //if((fixZ[s] == 1) & neut) continue;
            if((fixZ[s] == 1) && (i % 3 == 1)) continue;
            
            tmp_diag  = exp(eigenvals[i]*distances2[s]);
            mat x =eigenvecs[i];
            x.each_col()%=tmp_diag;
            log_cache_TM_tmp[i][s] = log(eigeninvs[i] * x);
        }
        
    }
}


vector<int> BPP_C::Update_ZB()
{
    vector<int> changedZ;
    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s =*it;
        //if(missing[s] ) continue; //
        log_prob_back[s] = log_emission[s];
        
        //            if(fixZ[s] ==1)
        //            {
        //                log_prob_back[s][2] = -INFINITY; //log_prob_back[s][0] = -INFINITY;
        //                log_prob_back[s][5] = -INFINITY;
        //            }
    }
    
     log_prob_back[root].fill(0);
    
    //message passing from bottom to top
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s = *it;
        int p = parent2[s];
        
        log_prob_back[p].subvec(0,2) += BPP::log_multi(log_TM_Int[s], log_prob_back[s]);  //matrix * vector
        log_prob_back[p].subvec(3,5) += BPP::log_multi(log_TM_Int[s], log_prob_back[s]);
        
        //cout << log_prob_back[p] <<endl;
    }
    
    //update Z and B from top to bottom Z[N] = 0, B[N] = 0
    int old;
    vec log_trans_p(6), trans_p(6), trans_B(2);
    
    for(vector<int>::iterator it =nodes.end()-2;it>=nodes.begin();it--)  //exclude root, doesn't include missing nodes
    {
        
        int s = *it;
        int p = parent2[s];
        
        //if(fixZ[s] == 1) continue;
        //if(fixZ[s] == 2 && (missing[s] && s<S)) continue;
        
        old = Z[s]; //Z[s] * 2 + isGB[s];
        if(Z[p]!=2)
        {
            //log_trans_p.subvec(0, 2) =log_TM_Int[s].unsafe_col(Z[p]) + log_prob_back[s].subvec(0, 2);
            log_trans_p =log_TM_Int[s].unsafe_col(Z[p]) + log_prob_back[s];
            trans_p = BPP::log_sample(log_trans_p);
            
            //cout << exp(log_y).t()<<endl;
            //cout << trans_p.t()<<endl;
            
            unsigned int n[6];
            gsl_ran_multinomial(RNG, 6, 1,trans_p.memptr(),n);
            int nn;
            for(nn=0; nn<6;nn++)
            {
                if(n[nn]>0) break;
            }
            
            Z[s] =nn % 3;
            if(!missing[s]) isGB[s] = (nn > 2);
            
            if(nn>5) {
                cout<< "Sample Z error" <<endl;
            }
            
        }
        else{
            Z[s] =2;
            trans_B.at(0) =log_prob_back[s].at(2); trans_B.at(1) =log_prob_back[s].at(5);
            trans_B = BPP::log_sample(trans_B);
            unsigned int n[2];
            gsl_ran_multinomial(RNG, 2, 1,trans_B.memptr(),n);
            int nn;
            for(nn=0; nn<2;nn++)
            {
                if(n[nn]>0) break;
            }
            if(!missing[s]) isGB[s] = (nn > 0);
            
            if(nn>1) {
                cout<< "Sample B error" <<endl;
            }
        }
        
        //if(old != Z[s] * 2 + isGB[s]) changedZ.push_back(s);
        if(old != Z[s]) changedZ.push_back(s);// only record Z[s]
        // cout << log_prob_back[s] <<endl;
        // cout<<trans_p<<endl;
        
    }
    return(changedZ);
}


void BPP_C::integrate_ZB(BPP &bpp, bool gb, bool neut, vector<mat> & Qs, vector<vector<double>> & log_emission_tmp, vector<vector<mat>> & log_cache_TM_tmp)
{
    
    update_cache_TM(Qs, log_cache_TM_tmp); // attention of Qs!!: "if((fixZ[s] == 1) && (i % 3 == 1))"s
    getEmission_update(neut, gb, log_emission_tmp, log_cache_TM_tmp, bpp, -1); //neut, gb: no use

    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s =*it;
        //if(missing[s] ) continue;
        log_prob_back[s] = log_emission_tmp[s];
        
        //            if(fixZ[s] ==1)
        //            {
        //                log_prob_back[s][2] = -INFINITY; //log_prob_back[s][0] = -INFINITY;
        //                log_prob_back[s][5] = -INFINITY;
        //            }
    }
    log_prob_back[root].fill(0);
    //message passing from bottom to top
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s = *it;
        int p = parent2[s];
        
        log_prob_back[p].subvec(0,2) += BPP::log_multi(log_TM_Int[s], log_prob_back[s] );  //matrix * vector
        log_prob_back[p].subvec(3,5) += BPP::log_multi(log_TM_Int[s], log_prob_back[s]);
        //cout << log_prob_back[p] <<endl;
    }
}


//// updating Z, return a vector of Z is changed, from top to bottom, for all the g
//vector<int>  BPP_C::Update_Z(int num_base)
//{
//
//    vector<int> changedZ;
//    
//    
//    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++) //N-1
//    {
//        int s = *it;
//        int p = parent2[s];
//    
//        double emission =0, emission2 = 0, emission3 =0 ;
//        if(missing[s])
//        {
//            //confined log_prob_back in initMCMC for s<S; otherwise log_emission ==0
//            if(s>=S) log_prob_back[s].fill(0);
//            continue;
//        }
//       
//        //sum over all g
//        for(int g=0;g <GG;g++)
//        {
//            if(Tg[g][s] == num_base) continue;
//            
//            if(Tg[g][s] ==-1){
//                
//                emission += ambiguousS_null[g][s][Tg[g][p]];
//                emission2 += ambiguousS_cons[g][s][Tg[g][p]];
//                emission3 += ambiguousS_neut[g][s][Tg[g][p]];
//                
//            }else{
//                emission += log_cache_TM_null[s].at(Tg[g][s],Tg[g][p]);
//                emission2 += log_cache_TM_cons[s].at(Tg[g][s],Tg[g][p]);
//                emission3 += log_cache_TM_neut[s].at(Tg[g][s],Tg[g][p]);
//            }
//        }
//        log_prob_back[s][0] =  emission; log_prob_back[s][1] =  emission2; log_prob_back[s][2] =  emission3;
//
//    }
//    
//
//    //message passing from bottom to top
//    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
//    {
//        int s = *it;
//        int p = parent2[s];
//       
//        log_prob_back[p] += BPP::log_multi2(log_TM_Int[s], log_prob_back[s]);  //matrix * vector
//        //cout << log_prob_back[p] <<endl;
//    }
//
//
//    //update Z from top to bottom Z[N] = 0
//    vec log_trans_p, trans_p;
//    for(vector<int>::iterator it = nodes.end() -2; it >= nodes.begin();it--) //N-2
//    {
//
//        int s = *it; //also sample missing Z!!
//        int p = parent2[s];
//        
//        int old = Z[s];
//        
//        if(Z[p]!=2)
//        {
//            log_trans_p = log_TM_Int[s].col(Z[p]);
//       
//
//            log_trans_p =log_trans_p + log_prob_back[s];
//            trans_p = BPP::log_sample(log_trans_p);
//
//            //cout << exp(log_y).t()<<endl;
//            //cout << trans_p.t()<<endl;
//
//            unsigned int n[3];
//           
//            gsl_ran_multinomial(RNG, 3, 1,trans_p.memptr(),n);
//            int nn;
//            for(nn=0; nn<3;nn++)
//            {
//                if(n[nn]>0) break;
//            }
//           
//
//            Z[s] = nn;
//
//            if(nn>2) {
//                throw runtime_error("Sample Z error");
//            }
//        }
//        else{
//            Z[s] = 2;
//        }
//
//        if(old != Z[s]) changedZ.push_back(s);
//
//       // cout << log_prob_back[s] <<endl;
//       // cout<<trans_p<<endl;
//
//    }
//    return(changedZ);
//
//}


vector<int>  BPP_C::Move_Z(int & propConf, int & revConf, int & changeZ){
    
    // propose new Z
    vector<int> propZ = vector<int>(N, 0);
    vector<int> increaseZ;
    vector<int> decreaseZ;
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it ++ )
    {
        int s = *it;
        int p = parent2[s];
        
        propZ[s] = Z[s] ;
        if(Z[s] == 2 && Z[p] == 1)
        {
            decreaseZ.push_back(s);
        }
        else if(Z[s] == 1 && Z[p] == 0)
        {
            int* cc = children2[s];
            if(cc[0]  == 2 || cc[1]  == 2) continue;
            decreaseZ.push_back(s);

        }
        else if(fixZ[s]!=1 && Z[s] == 1 && Z[p] == 1)
         {
             int* cc = children2[s];
             if(s<S || (cc[0]  == 2 && cc[1]  == 2)){
                 increaseZ.push_back(s);
             }
         }
        else if(fixZ[s]!=0 && Z[s] == 0) //= moveroot
        {
            int* cc = children2[s];
            if(cc[0]  == 1 && cc[1]  == 1){
                increaseZ.push_back(s);
            }
        }
    }
    
    int deN = decreaseZ.size();
    int inN = increaseZ.size();
    
    int totN = deN + inN;
    propConf = totN;
    
    int ind = rand() % (deN + inN);
    
    // compute reverse prob: propZ -> Z
    if(ind < deN)
    {
        changeZ = decreaseZ[ind];
        propZ[changeZ] -= 1;
        int p = parent2[changeZ];
        //int* cc = children2[s];
        if(Z[changeZ] == 1)
        {
            // s : 1 -> 0, p can't 0 -> 1
            if(fixZ[p]!=0 && Z[children2[p][0]] == 1 && Z[children2[p][1]] == 1 ) //p!= moveroot (and root)
            {
                totN -= 1;
            }
            // s : 1 -> 0, c can't 1 -> 2
            if(changeZ>=S)
            {
                for(int i  = 0; i < 1;i ++)
                {
                    int c = children2[changeZ][i];
                    if(fixZ[c] == 1) continue;
                    if(c<S || (Z[children2[c][0]]  == 2 && Z[children2[c][1]]  == 2))
                    {
                        totN--;
                    }
                    
                    if(c<S || (Z[children2[c][0]]  == 1 && Z[children2[c][1]]  == 1))
                    {
                        totN++;
                    }
                }
            }
            
        }
        else if(Z[changeZ] == 2)
        {
            // s : 2 -> 1, p can't 1 -> 2
            if( Z[parent2[p]] == 1 && Z[children2[p][0]]  == 2 && Z[children2[p][1]]  == 2 && fixZ[p]!=1){
                totN--;
            }
            else if(Z[parent2[p]] == 0 && propZ[children2[p][0]]  == 1 && propZ[children2[p][1]]  == 1) // s : 2 -> 1, p can 1 -> 0
            {
                totN++;
            }
            // s : 2 -> 1, p can be 2 -> 1
            if(changeZ >= S)
            {
                totN += 2;
            }
        }
    }else{
        changeZ = increaseZ[ind - deN];
        int p = parent2[changeZ];
        propZ[changeZ] += 1;
        if(Z[changeZ] == 1){
            // s : 1 -> 2, p can't 1 -> 0
            if(Z[parent2[p]] == 0 && Z[children2[p][0]] == 1 && Z[children2[p][1]] == 1)
            {
                totN--;
            }else if(Z[parent2[p]] == 1 && propZ[children2[p][0]] == 2 && propZ[children2[p][1]] == 2 && fixZ[p]!=1)
            {
                totN++;
            }
            // s : 1 -> 2, c can't 2 -> 1
            if(changeZ >= S)
            {
                totN -= 2;
            }
        }else if(Z[changeZ] == 0){
            // s : 0 -> 1, p can 0 -> 1
            if(fixZ[p]!=0 && Z[children2[p][0]] == 1 && Z[children2[p][1]] == 1 ) //p != moveroot
            {
                totN++;
            }
            if(changeZ >= S)
            {
                // s : 0 -> 1, c can't 1 -> 0
                // s : 0 -> 1, c can 1 -> 2
                for(int i  = 0; i < 1;i ++)
                {
                    int c = children2[changeZ][i];
                    if( c<S || (Z[children2[c][0]]  == 1 && Z[children2[c][1]]  == 1)) totN --;
                    if(fixZ[changeZ]!=1 &&( c<S || (Z[children2[c][0]]  == 2 && Z[children2[c][1]]  == 2))) totN++;
                }
            }


        }

    }
    
    revConf = totN;
    
    return(propZ);
    
}



//double BPP_C::Update_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut, vector<mat> & log_cache_TM_cons, vector<bool> & visited)
//{
//    double result =0;
//    // 1. sending the lambda msg from leaves bottom up through the network
//    for(vector<int>::iterator it = internal_nodes.begin(); it < internal_nodes.end(); it++)//int s=S; s<=root; s++)
//    {
//        int s = *it;
//        if(visited[s] || missing[s]) continue;
//        for(int g=0;g<GG;g++)
//        {
//
//            if(Tg[g][s]==num_base) continue;
//
//            int* p = children2[s];
//            //lambda[g][s].fill(0);
//
//            for(int cc=0;cc<2;cc++)
//            {
//                int chi = p[cc];
//                assert(chi != -1);
//                if(missing[chi] || Tg[g][chi]==num_base) continue;
//
//                if(Z[chi] ==2) lambda[g][s] +=  BPP::log_multi(log_cache_TM_neut[chi],lambda[g][chi]);
//                else if (Z[chi] ==0) lambda[g][s] +=  BPP::log_multi(log_cache_TM_null[chi],lambda[g][chi]);
//                else lambda[g][s] +=  BPP::log_multi(log_cache_TM_cons[chi],lambda[g][chi]);
//            }
//        }
//
//    }
//
//    // 2. processing the distribution of root species
//    for(int g=0;g<GG;g++)
//    {
//        if(!visited[root]) lambda[g][root] += log_pi;
//        result += BPP::log_exp_sum(lambda[g][root]);
//    }
//
//    for(vector<int>::iterator it = nodes.begin(); it < nodes.end(); it++)
//    {
//        if(missing[*it]) result += log_emission[*it][Z[*it]];// for missing
//    }
//
//    result += prior_Z_subtree(Z);
//    return(result);
//
//}


    




