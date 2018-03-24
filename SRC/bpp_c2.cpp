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


void BPP_C::Eval2(BPP&bpp, int resZ)  // VB lower bound
{
    
    //compare chib method with VB lower bound
    // compute posterior of Z
    Z = Max_Z; //consensus Z
    
    vector<int> posterior_Z = vector<int>(N,0);
    vector<unsigned int> uncertainZ;
    vector< vector<int> > configZ;
    vector< int > numConfigZ;
    // find uncertain nodes
    for(int s =0; s<N;s++)
    {
        for(std::size_t i = num_burn; i < trace_Z.size(); i++)
        {
            if(trace_Z[i][s] == Max_Z[s]) posterior_Z[s]++;
        }
        
        if(posterior_Z[s] <= 0.99 * num_mcmc)
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
    unsigned int numConf = 1;
    
    //double max_nrate = trace_n_rate[Max_m];
    //double max_crate = trace_c_rate[Max_m];
    
    vector<double> var_nrate = vector<double>(1,0);
    vector<double> mean_nrate = vector<double>(1,0);  // or mode?
    vector<double> var_crate = vector<double>(1,0);
    vector<double> mean_crate = vector<double>(1,0);
    vector<double> cov_rate = vector<double>(1,0);
    
    vector<double> sum_full_loglik = vector<double>(1,0);
    //double sum_full_loglik;
    //vector<double> max_full_loglik = vector<double>(1,-INFINITY);  //for chib
    
    //rate distribution of config MaxZ
    //vector<double> nrate;
    //vector<double> crate;
    //int checkConf = 0;
    for(int i = trace_Z.size() - 1; i >= num_burn ; i--)
    {
        unsigned int j = 0;
        for(; j < numConf;j++)
        {
            vector<unsigned int>::iterator it;
            for(it = uncertainZ.begin(); it <uncertainZ.end(); it++)
            {
                if(trace_Z[i][*it] != configZ[j][*it]){
                    
                    break;
                }
            }
            
            if(it == uncertainZ.end())
            {
               
                break;
            }
            
        }
        
        if(j == numConf)
        {
            configZ.push_back(trace_Z[i]);
            numConf++;
            numConfigZ.push_back(0);
            var_crate.push_back(0);
            var_nrate.push_back(0);
            mean_crate.push_back(0);
            mean_nrate.push_back(0);
            cov_rate.push_back(0);
            sum_full_loglik.push_back(0); //
            //max_full_loglik.push_back(-INFINITY); //

        }//else{
        
            //if(!ctnutils::IsVectorSame(trace_Z[i], configZ[j])) checkConf++;
        //}

        numConfigZ[j] ++;
        var_crate[j] += pow(trace_c_rate[i], 2);
        mean_crate[j] += trace_c_rate[i];
        if(resZ!=0){
            var_nrate[j] += pow(trace_n_rate[i], 2);  //no nrate for null model
            mean_nrate[j] += trace_n_rate[i];
            cov_rate[j] += trace_c_rate[i] * trace_n_rate[i];
        }
        sum_full_loglik[j] += trace_full_loglik[i]; //[j]
        //if(trace_full_loglik[i] > max_full_loglik[j])  max_full_loglik[j] = trace_full_loglik[i];
        
    }
    
    // H(q(Z))
    double VB =0; // sum_full_loglik/num_mcmc;
    double logdet;
    
    int effective_sum = 0;
    
    // get log(P(X, Z)),logdet(Sigma) for each configuration
    for(unsigned int i = 0; i<numConf;i++)
    {
        double w = (double)numConfigZ[i]/num_mcmc;
        if(w < 0.1) continue; // always keep MaxZ (if only one config, same as Chib)
        
        effective_sum += numConfigZ[i];
        
        //VB -= w * log(w);
        VB -= numConfigZ[i] * log(numConfigZ[i]);
        
        var_crate[i] -= pow(mean_crate[i],2)/numConfigZ[i] ;
        var_nrate[i] -= pow(mean_nrate[i],2)/numConfigZ[i] ;
        cov_rate[i] -= mean_nrate[i] * mean_crate[i]/numConfigZ[i] ;
        
        if(resZ==0)  // no or few Z==2, var_nrate[i] < 1e-6
        {
            logdet = log(var_crate[i]) - log(numConfigZ[i]);
        }else{
            logdet = log(var_crate[i] * var_nrate[i] - cov_rate[i] * cov_rate[i]) - 2*log(numConfigZ[i]);
        }
        
        if(verbose) cout << CC << ": weight: " << w << " logdet: " << logdet <<" "<<sum_full_loglik[i]/numConfigZ[i] + 0.5*logdet + log(2*M_PI) + 1 - log(w) << endl;

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
    if(effective_sum == 0 || logdet == NAN)
    {
        double varn =0 , varc = 0, cov = 0, meann = 0, meanc = 0;
        VB = 0;
        for(unsigned int i = 0; i < numConf;i++)
        {
            varn += var_nrate[i];
            varc += var_crate[i];
            cov += cov_rate[i];
            meann += mean_nrate[i];
            meanc += mean_crate[i];
            VB += sum_full_loglik[i];
            VB -= numConfigZ[i] * log(numConfigZ[i]);
        }
        if(resZ==0)  // no or few Z==2, var_nrate[i] < 1e-6
        {
            varc -= pow(meanc,2)/num_mcmc;
            logdet = log(varc) - log(num_mcmc);
        }else{
            varn -= pow(meann,2)/num_mcmc;
            varc -= pow(meanc,2)/num_mcmc;
            cov -= meann * meanc/num_mcmc ;
            logdet = log(varc * varn - cov * cov) - 2*log(num_mcmc);
            
        }
        
        VB += 0.5*logdet*num_mcmc;
        effective_sum = num_mcmc;
    }
    
    //chib method
    if(resZ==1) // full
    {
        if(effective_sum==0)
        {
          bpp.log_liks_sgl[CC] = MaxLoglik;
        }else{
          bpp.log_liks_sgl[CC] = VB/effective_sum + log(effective_sum) + log(2*M_PI) + 1;
        }
        if(verbose) cout << "P(X,Z*,r*): " << MaxLoglik << " P(X|null): " <<  bpp.log_liks_null[CC] <<" P(X|resZ): " << bpp.log_liks_resZ[CC] << " P(X): " << bpp.log_liks_sgl[CC] << endl;
        //cout << "P(X,Z*,r*): " << MaxLoglik << " P(X): " << bpp.log_liks_sgl[CC] << endl;

        if(bpp.log_liks_sgl[CC] == NAN) failure = 1;
    }
    else if(resZ==0)
    {
        bpp.log_liks_null[CC] = VB/effective_sum + log(effective_sum) + 0.5*(log(2*M_PI) + 1);
        if(bpp.log_liks_null[CC] == NAN) failure = 1;

    }else{
        bpp.log_liks_resZ[CC] = VB/effective_sum + log(effective_sum) + log(2*M_PI) + 1;
        if(bpp.log_liks_resZ[CC] == NAN) failure = 1;
    }
    
    
    
}

void BPP_C::getEmission(int num_base)
{
    
    for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++) //N-1
    {
        int s = *it;
        int p = parent2[s];
        
        if(missing[s]) continue; //log_emission is zero, no!! already initialized at initMCMC!
        
        //log_TM = log_cache_TM_null[s];
        double emission =0, emission2 = 0, emission3 =0 ;
        //sum over all g
        for(int g=0;g <GG;g++)
        {
            if(Tg[g][s] == num_base) continue;
            
            if(Tg[g][s] ==-1){
                emission += ambiguousS_null[g][s][Tg[g][p]];
                //emission2 += ambiguousS_cons[g][s][Tg[g][p]];
                //emission3 += ambiguousS_neut[g][s][Tg[g][p]];
                emission2 += BPP::log_exp_sum(log_cache_TM_cons[s].col(Tg[g][p]) + lambda[g][s]);
                emission3 += BPP::log_exp_sum(log_cache_TM_neut[s].col(Tg[g][p]) + lambda[g][s]);
                
            }else{
                emission += log_cache_TM_null[s].at(Tg[g][s],Tg[g][p]);
                emission2 += log_cache_TM_cons[s].at(Tg[g][s],Tg[g][p]);
                emission3 += log_cache_TM_neut[s].at(Tg[g][s],Tg[g][p]);
            }
        }
        log_emission[s][0] =  emission; log_emission[s][1] =  emission2; log_emission[s][2] =  emission3;
        
    }
    
    
}

double BPP_C::log_f_Xz(vec log_pi, int num_base, vector<int>& Z, vector<mat> & log_cache_TM_neut,  vector<mat> & log_cache_TM_cons)
{
    double result =0;
    for(int g=0;g<GG;g++)
    {
        // 1. sending the lambda msg from leaves bottom up through the network
        for(vector<int>::iterator it = internal_nodes.begin(); it < internal_nodes.end(); it++)//int s=S; s<=root; s++)
        {
            int s = *it;
            if(missing[s] || Tg[g][s]==num_base) continue;
            
            int* p = children2[s];
            lambda[g][s].fill(0);
            
            for(int cc=0;cc<2;cc++)
            {
                int chi = p[cc];
                assert(chi != -1);
                if(missing[chi] || Tg[g][chi]==num_base) continue;
                
                if(Z[chi] ==2) lambda[g][s] +=  BPP::log_multi(log_cache_TM_neut[chi],lambda[g][chi]);
                else if (Z[chi] ==0) lambda[g][s] +=  BPP::log_multi(log_cache_TM_null[chi],lambda[g][chi]);
                else lambda[g][s] +=  BPP::log_multi(log_cache_TM_cons[chi],lambda[g][chi]);
            }
            
        }
        
        // 2. processing the distribution of root species
        
        lambda[g][root] += log_pi;
        result +=BPP::log_exp_sum(lambda[g][root]);
    }
    return(result);
}


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
        
        log_prob_back[p] += BPP::log_multi2(log_TM_Int[s], log_prob_back[s]);  //matrix * vector
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
    
    #pragma omp critical
    {
      if(iter ==0)
      {
        out_lik.open(outpath_lik.c_str());
        out_lik << "loglik\trate_n\trate_c\t";
        for(int s =0 ;s<N;s++){  // header: species name
            out_lik<<bpp.nodes_names[s] << "\t";
        }
        out_lik <<endl;
      }else{
        out_lik.open(outpath_lik.c_str(), ios::app);
      }
        
      for(std::size_t i=0; i< trace_loglik.size();i++)
        {
            out_lik<<trace_loglik[i]<<"\t"<<trace_n_rate[i] << "\t"<<trace_c_rate[i]<<"\t";
            for(int s=0; s<N;s++)
                out_lik<<trace_Z[i][s]<<"\t";
            out_lik <<endl;
        }
    
    out_lik.close();
    }
}


void BPP_C::Output_init(string output_path,string output_path2, BPP &bpp, ofstream & out_Z, int resZ)
{
    
    std::sort(trace_n_rate.begin(), trace_n_rate.end());
    double n_rate =  trace_n_rate[trace_n_rate.size()/2];
    
    std::sort(trace_c_rate.begin(), trace_c_rate.end());
    double c_rate =  trace_c_rate[trace_c_rate.size()/2];
    
    vector<vector<int>> countZ = vector<vector<int>> (N,vector<int>(4,0));
    for(int s=0; s<N;s++)
    {
        
        for(std::size_t i = num_burn; i< trace_loglik.size();i++)
        {
            
            countZ[s][trace_Z[i][s]+1]++;
            
        }
        
        if(missing[s])
        {
            countZ[s][0] = num_mcmc;  //set missing s = 1, though Z[s] can be 0/1/2; only missing in upper Z[s] = -1
        }
    }
    
    #pragma omp critical
    {
        out_Z << CC << "\t" << n_rate << "\t" << c_rate;
        for(int s=0; s<N;s++){
            //out_Z <<"\t"<<countZ[s][0];
            for(int k=0;k<4;k++)
                out_Z <<"\t"<<(double)countZ[s][k]/num_mcmc; //?-1
        }
        
        out_Z <<endl;
    }
}









