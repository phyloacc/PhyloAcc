//
//  bpp_c2.cpp
//  PhyloAcc_init2
//
//  Created by hzr on 4/20/16.
//  Copyright © 2016 hzr. All rights reserved.
//

#include "bpp_c.hpp"
#include "bpp.hpp"
#include <iomanip>


set<int> BPP_C:: getShifts(vector<bool> & Z)
{
    set<int> res;
    for(vector<int>::iterator it = nodes.begin(); it <nodes.end(); it++)
    {
        if(missing[*it]) continue;
        if(Z[*it])
        {
            res.insert(*it);
        }
    }
    
    return(res);
}

void BPP_C::Eval2(BPP&bpp, int resZ)  // VB lower bound
{

    //compare chib method with VB lower bound
    // compute posterior of Z
    Z = Max_Z; //consensus Z

    vector<int> posterior_Z = vector<int>(N,0);
    vector<unsigned int> uncertainZ;
    vector< vector<int> > configZ;
    vector<vector<set<int>>> configZB; // record the isGB nodes for each config of Z
    vector< int > numConfigZ;
    vector< vector<int >> numConfigZB;
    map<set<int>, int> numConfigB; // overall config of B
    // find uncertain nodes
    for(int s =0; s<N;s++)
    {
        for(std::size_t i = num_burn; i < trace_Z.size(); i++)
        {
            if(trace_Z[i][s] == Max_Z[s]) posterior_Z[s]++;
        }

        if(posterior_Z[s] < num_mcmc)
        {
            uncertainZ.push_back(s);
        }
    }

//    cout << "uncertainZ: ";
//    ctnutils::DispVector(uncertainZ);
//    cout << endl;

    //find configurations of Z
    configZ.push_back(Max_Z);
    numConfigZ.push_back(0);
    numConfigZB.push_back(vector<int>());
    configZB.push_back(vector<set<int>>());
    unsigned int numConf = 1;

    int d = 2 + (resZ!=0); // dimension of MVN
    
    vector<vec> means = vector<vec>(1,zeros<vec>(d));  // or mode?, gB, crate and nrate (no for resZ == 0)
    vector<mat> var = vector<mat>(1,zeros<mat>(d, d)); // covariates matrix
    vector<double> sum_full_loglik = vector<double>(1, 0.0);//vector<vector<double>>(1,vector<double>(1, 0.0));
   
    for(int i = trace_Z.size() - 1; i >= num_burn ; i--)
    {
        unsigned int j = 0;
        for(; j < numConf;j++)
        {
            vector<unsigned int>::iterator it;
            for(it = uncertainZ.begin(); it <uncertainZ.end(); it++)
            {
                if(trace_Z[i][*it] != configZ[j][*it]) break;

            }
            if(it == uncertainZ.end())  break;
        }

        if(j == numConf)
        {
            configZ.push_back(trace_Z[i]);
            numConf++;
            numConfigZ.push_back(0);
            numConfigZB.push_back(vector<int>());
            configZB.push_back(vector<set<int>>());
            means.push_back(zeros<vec>(d));
            var.push_back(zeros<mat>(d,d));
            sum_full_loglik.push_back(0); //vector<double>()
        }
        numConfigZ[j] ++;
        
        
        set<int> tempSet = getShifts(trace_GB[i]) ;
        
        map<set<int>,int>::iterator it = numConfigB.find(tempSet);  // update numConfigB
        if(it!=numConfigB.end())
        {
             ++ it->second;
        }else{
            numConfigB.insert(pair<set<int>,int> (tempSet,1));
        }
        
        unsigned int k = 0;
        for(; k < numConfigZB[j].size();k ++)
        {
            if(tempSet == configZB[j][k]) break;
        }
        
        if(k == numConfigZB[j].size())
        {
            configZB[j].push_back(tempSet);
            numConfigZB[j].push_back(0);
            //sum_full_loglik[j].push_back(0.0);
        }
        numConfigZB[j][k] ++ ;
        
        var[j].at(1,1) += pow(trace_c_rate[i], 2);
        means[j].at(1) += trace_c_rate[i];

        var[j].at(0,0) += pow(trace_n_GB[i], 2);
        means[j].at(0) += trace_n_GB[i];
        
        var[j].at(0,1) += trace_c_rate[i] * trace_n_GB[i];

        if(resZ!=0){
            var[j].at(2,2) += pow(trace_n_rate[i], 2);
            means[j].at(2) += trace_n_rate[i];
            var[j].at(1,2) += trace_c_rate[i] * trace_n_rate[i];
            var[j].at(0,2) += trace_n_GB[i] * trace_n_rate[i];
        }
        sum_full_loglik[j] += trace_full_loglik[i]; //[j]
        //if(trace_full_loglik[i] > max_full_loglik[j])  max_full_loglik[j] = trace_full_loglik[i];

    }

    // H(q(Z))
    double VB =0; // sum_full_loglik/num_mcmc;
    double logdet;
    double sign;
    int effective_sum = 0;
    
    // get log(P(X, Z)),logdet(Sigma) for each configuration
    for(unsigned int i = 0; i<numConf;i++)
    {
        double w = (double)numConfigZ[i]/num_mcmc;
        if(w < 0.05) continue; //0.1 always keep MaxZ (if only one config, same as Chib)

        mat temp_var = symmatu(var[i]);  // get the symmetric part of covariance
        temp_var -= means[i] * means[i].t()/numConfigZ[i];
        temp_var /= (numConfigZ[i] - 1);

        if(numConfigZ[i] > 10)
        {
            log_det( logdet, sign, temp_var );
        }

        if(numConfigZ[i] <= 10 || sign < 0 || logdet == NAN)
        {
            logdet = sum(log(temp_var.diag())); // q = indep norm
        }

        if(logdet == NAN || logdet == -INFINITY) continue;

        effective_sum += numConfigZ[i];
        //VB -= w * log(w);
        //VB -= numConfigZ[i] * log(numConfigZ[i]);

//        var[i].at(1,1) -= pow(means[i].at(1),2)/numConfigZ[i] ;
//        var[i].at(2,2) -= pow(means[i].at(2),2)/numConfigZ[i] ;
//        cov_rate[i] -= mean_nrate[i] * mean_crate[i]/numConfigZ[i]
//        if(resZ==0)  // no or few Z==2, var_nrate[i] < 1e-6
//        {
//            logdet = log(var_crate[i]) - log(numConfigZ[i]);
//        }else{
//            logdet = log(var_crate[i] * var_nrate[i] - cov_rate[i] * cov_rate[i]) - 2*log(numConfigZ[i]);
//        }
        
        // vector<int>::iterator max_ind = max_element(numConfigZB[i].begin(), numConfigZB[i].end()); // get the mode of B with Z
        // int k = distance(numConfigZB[i].begin(), max_ind);
        // get q(B|X, Z)
        double tempB = 0;
        for(vector<int>::iterator it = numConfigZB[i].begin(); it < numConfigZB[i].end(); it++)
        {
            tempB += (double)(*it) * log(*it);
        }
        VB -= tempB;
        
        //VB -= (double)numConfigZ[i] * log(numConfigZ[i]); //canceled!
        //cout << tempB << endl;
        //cout << temp_var << endl;
        if(verbose) cout << CC << ": weight: " << w << " logdet: " << logdet <<" "<<sum_full_loglik[i]/numConfigZ[i] + 0.5*logdet + (double)d/2 * log(2*M_PI) + (double)d/2*(numConfigZ[i] - 1)/numConfigZ[i] - tempB/numConfigZ[i] + log(num_mcmc) << endl;  // a little different! num_mcmc vs. effective_sum

        VB += (double)(numConfigZ[i] - 1)/2 * d;
        VB += sum_full_loglik[i];
        // compare with chib;
        //cout << "P(x, r* Z): " <<max_full_loglik[i]<<"  avg P(x, r* Z) + 1: " << sum_full_loglik[i]/numConfigZ[i] + 1 << endl;
        //compare P(X|Z)/P(Z|X) for each configuration
        VB += 0.5*logdet*numConfigZ[i];
    }
    //cout << endl;
    if(verbose) cout << CC << ": Total Num of different confZ: " << numConf << " Number of MaxZ:" << numConfigZ[0] << endl; //" Conf different: " << checkConf <<

    // if effective_sum == 0 (no major configuration of Z) or logdet == NAN (varince stick at one value);
    // then assume q(r, Z) = q(r) q(Z), r and Z are indepentdent; use all rate to get q(r) and all configuration to get Z
    double tempB = 0;
    if(effective_sum == 0)
    {
        vec mean0 = zeros<vec>(d);
        mat var0 = zeros<mat>(d, d);
        VB = 0;
        for(unsigned int i = 0; i < numConf;i++)
        {
            mean0 += means[i];
            var0 += var[i];
            VB += sum_full_loglik[i];
            //VB -= numConfigZ[i] * log(numConfigZ[i]);
            
            for(vector<int>::iterator it = numConfigZB[i].begin(); it < numConfigZB[i].end(); it++)
            {
                tempB += (double)(*it) * log(*it);
            }
           
            
        }
        VB -= tempB;
        //cout << tempB << endl;
        
//        for(map<set<int>, int> :: iterator it = numConfigB.begin(); it != numConfigB.end(); it++)
//        {
//            VB -= it->second * log(it->second);
//        }
        
        var0 = symmatu(var0);  // get the symmetric part of covariance
        var0 -= mean0 * mean0.t()/num_mcmc;
        var0 /= (num_mcmc - 1);
        
        
       // cout << var0 << endl;
        
        log_det( logdet, sign, var0);
        
        if(sign < 0 || logdet == NAN)
        {
            logdet = sum(log(var0.diag())); // q = indep norm
        }
        
        if(logdet != NAN && logdet != -INFINITY)
        {
            VB += 0.5*logdet*num_mcmc;
            //VB += log(num_mcmc) * num_mcmc;
            VB += (double)(num_mcmc - 1)/2 * d;
            effective_sum = num_mcmc;
        }
    }

    if(resZ==1) // full
    {
        if(effective_sum==0)
        {
          bpp.log_liks_sgl[CC] = MaxLoglik;
        }else{
            bpp.log_liks_sgl[CC] = VB/effective_sum + log(effective_sum) + (double)d/2*log(2*M_PI); //+ 1;
        }
        if(verbose) cout << "P(X,Z*,r*): " << MaxLoglik << " P(X|null): " <<  bpp.log_liks_null[CC] <<" P(X|resZ): " << bpp.log_liks_resZ[CC] << " P(X): " << bpp.log_liks_sgl[CC] << endl;
        //cout << "P(X,Z*,r*): " << MaxLoglik << " P(X): " << bpp.log_liks_sgl[CC] << endl;
    }
    else if(resZ==0)
    {
        bpp.log_liks_null[CC] = VB/effective_sum + log(effective_sum) + (double)d/2*log(2*M_PI);
    }else{
        bpp.log_liks_resZ[CC] = VB/effective_sum + log(effective_sum) + (double)d/2*log(2*M_PI) ;
    }



}

void BPP_C::getEmission_update(bool neut, bool gb, vector<vector<double>> & log_emission_tmp,  vector<vector<mat >>& log_cache_TM_tmp, BPP &bpp, int restore)
{
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++) //N-1
    {
        int s = *it;
        if(missing[s]) continue; //already initialized at initMCMC!
        
        mat transition;
        
        vector<int> inds; // which coordinates in log_emission to be updated
//        if(gb)
//        {
//            inds.push_back(3);   inds.push_back(4);   inds.push_back(5);
//        }else if(!neut)
//        {
//            inds.push_back(1);   inds.push_back(4);
//        }else{
//            inds.push_back(2);   inds.push_back(5);
//        }
        for(int i =1; i < 6; i ++) inds.push_back(i);
        
        for(int i =0 ;i <inds.size(); i++)
        {
            int j  = inds[i];
            if(fixZ[s] == 1 && (j % 3 ==2)) continue; //-infinity
            if(restore == -1)
            {
                log_emission_tmp[s][j]  = accu(log_cache_TM_tmp[i][s] % posterior_Tg[s]);
                
                if(j < 3) log_emission_tmp[s][j] += log( 1- bpp.nGBprior_c) ;
                else  log_emission_tmp[s][j] += log(bpp.nGBprior_c);
                
                if(j % 3 == 1) log_emission_tmp[s][j] += log(1- bpp.consToMis);
                else log_emission_tmp[s][j] += log(1- bpp.nconsToMis);
            }
            else if(restore == 0) // MH rejected
            {
                log_emission_tmp[s][j] = log_emission[s][j];
            }else // MH accepted
            {
                log_emission[s][j] = log_emission_tmp[s][j];
                switch (j) {
                    case 1:
                        log_cache_TM_cons[s] = log_cache_TM_tmp[i][s];
                        break;
                    case 2:
                        log_cache_TM_neut[s] = log_cache_TM_tmp[i][s];
                        break;
                    case 3:
                        log_cache_TM_null_B[s] = log_cache_TM_tmp[i][s];
                        break;
                    case 4:
                        log_cache_TM_cons_B[s] = log_cache_TM_tmp[i][s];
                        break;
                    case 5:
                        log_cache_TM_neut_B[s] = log_cache_TM_tmp[i][s];
                        break;
                    default:
                        break;
                }
            }
        }
        
    }
    
}

double BPP_C::initial_prob_back(BPP &bpp) // update log_prob_back within
{
    double elbo = 0;
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++) //N-1
    {
        int s = *it;
       // int p = parent2[s];
        
        if(missing[s]) continue; //log_emission is zero, no!! already initialized at initMCMC!
        
        mat transition;
        for(int i =0 ;i <6; i++)
        {
            
            if(fixZ[s] == 1 && (i % 3 ==2)) continue; //-infinity
            switch (i) {
                case 0:
                    transition = log_cache_TM_null[s]; break;
                case 3:
                    transition = log_cache_TM_null_B[s]; break;
                case 1:
                    transition = log_cache_TM_cons[s];
                    break;
                case 4:
                    transition = log_cache_TM_cons_B[s];
                    break;
                case 2:
                    transition = log_cache_TM_neut[s];
                    break;
                case 5:
                    transition = log_cache_TM_neut_B[s];
                    break;
            }
            log_emission[s][i]  = accu(transition % posterior_Tg[s]);
            
            if( i== (Z[s] + 3 *isGB[s])) elbo += log_emission[s][i] ;
            
            if(i < 3) log_emission[s][i] += log( 1- bpp.nGBprior_c) ;
            else  log_emission[s][i] += log(bpp.nGBprior_c);
            
            if(i % 3 == 1) log_emission[s][i] += log(1- bpp.consToMis);
            else log_emission[s][i] += log(1- bpp.nconsToMis);
        }
        
    }
    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s =*it;
        //if(missing[s] ) continue; 
        log_prob_back[s] = log_emission[s];
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
    return(elbo);
}

//double BPP_C::log_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut,  vector<mat> & log_cache_TM_cons)
//{
//    double result =0;
//    for(int g=0;g<GG;g++)
//    {
//        // 1. sending the lambda msg from leaves bottom up through the network
//        for(vector<int>::iterator it = internal_nodes.begin(); it < internal_nodes.end(); it++)//int s=S; s<=root; s++)
//        {
//            int s = *it;
//            if(missing[s] || Tg[g][s]==num_base) continue;
//
//            int* p = children2[s];
//            lambda[g][s].fill(0);
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
//
//        }
//
//        // 2. processing the distribution of root species
//
//        lambda[g][root] += log_pi;
//        result +=BPP::log_exp_sum(lambda[g][root]);
//    }
//    return(result);
//}


void BPP_C::log_f_Z(vector<int>& Z, vector<mat> & log_Int, double & MH_ratio_g, double & MH_ratio_l)
{
    MH_ratio_g = 0;
    MH_ratio_l = 0;
    
    
    for(vector<int>::iterator it = nodes.begin(); it <nodes.end() - 1; it++)
    {
        int s = *it;
        if(Z[parent2[s]] == 0)
        {
            if(log_TM_Int[s](1,0) == 0 || log_TM_Int[s](0,0) == 0) continue;
            MH_ratio_g +=  log_TM_Int[s](1,0) * Z[s] + log_TM_Int[s](0,0) * (1 - Z[s]);
        }else if(Z[parent2[s]] == 1){
            if(log_TM_Int[s](2,1) == 0 || log_TM_Int[s](1,1) == 0) continue;
            MH_ratio_l +=  log_TM_Int[s](2,1) * (Z[s] - 1) + log_TM_Int[s](1,1) * (2 - Z[s]);
        }
        
        
    }
    
    
}


double BPP_C::prior_Z(vector<int> & tmpZ)  // no constraint
{
    double result =0;
    
    //update Z from top to bottom Z[N] = 0
    for(vector<int>::iterator it =nodes.end()-2;it>=nodes.begin();it--)  //
    {
        
        int s = *it;
        int p = parent2[s];
        result += log_TM_Int[s](tmpZ[s], tmpZ[p]);
   
    }
    return(result);
}

double BPP_C::prior_Z_subtree(vector<int> & tmpZ)  //whether cal prob_back again, default cal=0
{
    double result =0;
    //vector<vec> log_prob_back_temp = vector<vec> (N, zeros<vec>(3));
    
    // restore prior of Z (no X, no missing only constraint!)
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s = *it;
        //if(missing[s] && s<S) continue; //don't update log_prob_back for missing nodes
        log_prob_back[s].fill(0);
        if(fixZ[s] ==1)
        {
            log_prob_back[s][2] = -INFINITY; //log_prob_back[s][0] = -INFINITY;
        }
        //        else if(fixZ[s] ==2)  // no use!
        //        {
        //            log_prob_back_temp[s][0] = -INFINITY;  //log_prob_back[s][1] = -INFINITY;
        //        }
    }
    
    
    //message passing from bottom to top
    mat tmp = zeros<mat> (3,2);
    
    //message passing from bottom to top
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
    {
        int s = *it;
        int p = parent2[s];
        
        log_prob_back[p] += BPP::log_multi(log_TM_Int[s], log_prob_back[s]);  //matrix * vector
        //cout << log_prob_back[p] <<endl;
    }
    
    
    
    //update Z from top to bottom Z[N] = 0
    vec log_trans_p, trans_p;
    for(vector<int>::iterator it =nodes.end()-2;it>=nodes.begin();it--)  //
    {
        
        int s = *it;
        int p = parent2[s];
        
        if(tmpZ[p]!=2)
        {
            
            log_trans_p = log_TM_Int[s].unsafe_col(tmpZ[p]) + log_prob_back[s];
            result += log_trans_p[tmpZ[s]] - BPP::log_exp_sum(log_trans_p);
            
            //trans_p = BPP::log_sample(log_trans_p);
            
            
        }else if(tmpZ[s]!=2){
            cout <<"prior_Z_subtree error, Z[s]!=2, p: " << p <<"s: "<<s <<endl;
        }
        
    }
    
    return(result);
    // cout << endl;
}




void BPP_C::Output_sampling(int iter, string output_path2, BPP &bpp, int resZ){
    
    string outpath_lik = output_path2 +"_mcmc_trace_" + to_string(resZ) + "_" + to_string(CC) +".txt";
    ofstream out_lik;
    out_lik.precision(8);
    
    if(iter ==0)
    {
        out_lik.open(outpath_lik.c_str());
        out_lik << "loglik\trate_n\trate_c\tGB\t";
        for(int s =0 ;s<N;s++){  // header: species name
            out_lik<<bpp.nodes_names[s] << "_B\t";
        }
        for(int s =0 ;s<N;s++){  // header: species name
            out_lik<<bpp.nodes_names[s] << "\t";
        }
        out_lik <<endl;
    }else{
        out_lik.open(outpath_lik.c_str(), ios::app);
    }
    #pragma omp critical
    {
        for(std::size_t i=0; i< trace_loglik.size();i++)
        {
            out_lik<<trace_full_loglik[i]<<"\t"<<trace_n_rate[i] << "\t"<<trace_c_rate[i]<<"\t"<<trace_n_GB[i] << "\t"; //<<trace_c_GB[i]<<"\t";
            for(int s=0; s<N;s++)
                out_lik<<trace_GB[i][s]<<"\t";
            for(int s=0; s<N;s++)
                out_lik<<trace_Z[i][s]<<"\t";
            out_lik <<endl;
        }
    }
    out_lik.close();
}


void BPP_C::Output_init(string output_path,string output_path2, BPP &bpp, ofstream & out_Z, int resZ)
{
    
    std::sort(trace_n_rate.begin() + num_burn, trace_n_rate.end());
    double n_rate =  trace_n_rate[trace_n_rate.size()/2];
    
    std::sort(trace_c_rate.begin() + num_burn, trace_c_rate.end());
    double c_rate =  trace_c_rate[trace_c_rate.size()/2];
    
    std::sort(trace_n_GB.begin() + num_burn, trace_n_GB.end());
    double gB =  trace_n_GB[trace_n_GB.size()/2];
    
    vector<vector<int>> countZ = vector<vector<int>> (N,vector<int>(4,0));
    vector<int> countB = vector<int> (N,0);
    for(int s=0; s<N;s++)
    {
        //vector<double> temp = vector<double>(num_mcmc, 0.0);
//        double sum = 0.0;
//        for(std::size_t i = num_burn; i< num_mcmc + num_burn;i++)
//        {
//            sum += trace_n_GB[i][s];
//        }
//        //std::sort(temp.begin() + num_burn, temp.end());
//        medianB[s] =  sum/num_mcmc; //temp[temp.size()/2];
//
        for(std::size_t i = num_burn; i< trace_loglik.size();i++)
        {
            
            countZ[s][trace_Z[i][s]+1]++;
            countB[s] += trace_GB[i][s];
        }
        
        if(missing[s])
        {
            countZ[s][0] = num_mcmc;  //set missing s = 1, though Z[s] can be 0/1/2; only missing in upper Z[s] = -1
        }
    }
    
    #pragma omp critical
    {
        out_Z << CC << "\t" << n_rate << "\t" << c_rate << "\t" << gB;
//        for(int s=0; s<N;s++){
//            out_Z <<"\t"<< medianB[s];
//        }
        for(int s=0; s<N;s++){
                out_Z <<"\t"<<(double)countB[s]/num_mcmc; //?-1
        }
        
        for(int s=0; s<N;s++){
            //out_Z <<"\t"<<countZ[s][0];
            for(int k=0;k<4;k++)
                out_Z <<"\t"<<(double)countZ[s][k]/num_mcmc; //?-1
        }
        
        
        
        out_Z <<endl;
    }
}









