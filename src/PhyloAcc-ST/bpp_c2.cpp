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
        sum_full_loglik[j] += trace_full_loglik[i]/num_mcmc; //[j]
        //if(trace_full_loglik[i] > max_full_loglik[j])  max_full_loglik[j] = trace_full_loglik[i];
        
    }
    
    // H(q(Z))
    double VB =0; // sum_full_loglik/num_mcmc;
    double logdet;
    
    int effective_sum = 0;
    
    //// get P(Z |C) for each configuration
    ////vector<double> priorZ = prior_Z_subtree(configZ, numConfigZ);
    // get log(P(X, Z)),logdet(Sigma) for each configuration
    for(unsigned int i = 0; i<numConf;i++)
    {
        double w = (double)numConfigZ[i]/num_mcmc;
        if(w < 0.02) continue; // always keep MaxZ (if only one config, same as Chib)
        
        // compute P(Z|C)
        mat nZ = zeros(3,2);
        for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
        {
            int p = parent2[*it];
            if(configZ[i][p] < 2)
            {
                nZ(configZ[i][*it],configZ[i][p]) += 1;
            }
        }
        
        double priorZ = gsl_sf_lnbeta(prior_g_a + nZ(1,0) + nZ(2,0), prior_g_b + nZ(0,0)) - gsl_sf_lnbeta(prior_g_a, prior_g_b);
        
        for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
        {
            if(fixZ[*it] == 1)
            {
                int p = parent2[*it];
                assert(p!=N);
                if(configZ[i][p] < 2)
                {
                    nZ(configZ[i][*it],configZ[i][p]) -= 1;
                }
            }
        }
        
        priorZ += gsl_sf_lnbeta(prior_l_a + nZ(2,1), prior_l_b + nZ(1,1)) - gsl_sf_lnbeta(prior_l_a, prior_l_b);
        if(prior_l2_a > 0) priorZ += gsl_sf_lnbeta(prior_l2_a + nZ(2,0), prior_l2_b + nZ(1,0)) - gsl_sf_lnbeta(prior_l2_a, prior_l2_b);
        
        //root
        priorZ += prior_z.at(configZ[i][root]);
        
        sum_full_loglik[i] += numConfigZ[i] * priorZ /num_mcmc;
        
        if(resZ==0)  // no or few Z==2, var_nrate[i] < 1e-6
        {
            logdet = log(var_crate[i] - pow(mean_crate[i],2)/numConfigZ[i]) - log(numConfigZ[i]);
        }else{
            double vc = var_crate[i] - pow(mean_crate[i],2)/numConfigZ[i];
            double vn = var_nrate[i] - pow(mean_nrate[i],2)/numConfigZ[i] ;
            double cov = cov_rate[i] - mean_nrate[i] * mean_crate[i]/numConfigZ[i] ;
            logdet = log(vc * vn - cov * cov) - 2*log(numConfigZ[i]);
        }
        
        if(::isnan(logdet)) continue;
        
        effective_sum += numConfigZ[i];
        
        //VB -= w * log(w);
        VB -=  w * log(numConfigZ[i]); //numConfigZ[i]
        
        VB += sum_full_loglik[i];
        
        // compare with chib;
        //cout << "P(x, r* Z): " <<max_full_loglik[i]<<"  avg P(x, r* Z) + 1: " << sum_full_loglik[i]/numConfigZ[i] + 1 << endl;
        //compare P(X|Z)/P(Z|X) for each configuration
        VB += 0.5*logdet*numConfigZ[i]/num_mcmc;
        
        if(verbose) cout << CC << ": weight: " << w << " logdet: " << logdet <<" "<<sum_full_loglik[i]/numConfigZ[i] * num_mcmc + 0.5*logdet + log(2*M_PI) + 1 - log(w) << endl;
        // the difference from final result is w is divided by num_mcmc not effective sum, so w is smaller, loglik is larger
    }
    //cout << endl;
    if(verbose) cout << CC << ": Total Num of different confZ: " << numConf << " Number of MaxZ:" << numConfigZ[0] << endl << endl ; //" Conf different: " << checkConf <<
    
    // if effective_sum == 0 (no major configuration of Z) or logdet == NAN (varince stick at one value);
    // then assume q(r, Z) = q(r) q(Z), r and Z are indepentdent; use all rate to get q(r) and all configuration to get Z
    if(effective_sum == 0 ) //|| logdet == NAN
    {
        double varn =0 , varc = 0, cov = 0, meann = 0, meanc = 0;
        VB = 0;
        for(unsigned int i = 0; i < numConf;i++)
        {
            
            varc += var_crate[i];
            meanc += mean_crate[i];
            if(resZ!=0){
                meann += mean_nrate[i];
                varn += var_nrate[i];
                cov += cov_rate[i];
            }
           
        }
        
        double rn = meann /num_mcmc; // get r*
        double rc = meanc /num_mcmc;
        
        // recompute log_cache_TM at r*
        for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
        {
            int s = *it;
            if(s==bpp.moveroot) continue; //not update rate for outgroup! s==42||, conserved rate still changing
            
            mat tmp_diag  = exp(bpp.eigenval*distances2[s]*rc);
            mat x = bpp.eigenvec;
            x.each_col()%=tmp_diag;
            log_cache_TM_cons[s] = log(bpp.eigeninv * x);
            
            if(resZ!=0)
            {
                tmp_diag  = exp(bpp.eigenval*distances2[s]*rn);
                x = bpp.eigenvec;
                x.each_col()%=tmp_diag;
                
                log_cache_TM_neut[s] = log(bpp.eigeninv * x); //transpose Q
            }
        }
        
        
        for(unsigned int i = 0; i < numConf;i++)
        {
            // compute P(Z|C)
            mat nZ = zeros(3,2);
            for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
            {
                int p = parent2[*it];
                if(configZ[i][p] < 2)
                {
                    nZ(configZ[i][*it],configZ[i][p]) += 1;
                }
            }
            double priorZ = gsl_sf_lnbeta(prior_g_a + nZ(1,0) + nZ(2,0), prior_g_b + nZ(0,0)) - gsl_sf_lnbeta(prior_g_a, prior_g_b);
            
            for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
            {
                if(fixZ[*it] == 1)
                {
                    int p = parent2[*it];
                    assert(p!=N);
                    if(configZ[i][p] < 2)
                    {
                        nZ(configZ[i][*it],configZ[i][p]) -= 1;
                    }
                }
            }
            
            priorZ += gsl_sf_lnbeta(prior_l_a + nZ(2,1), prior_l_b + nZ(1,1)) - gsl_sf_lnbeta(prior_l_a, prior_l_b);
            if(prior_l2_a > 0) priorZ += gsl_sf_lnbeta(prior_l2_a + nZ(2,0), prior_l2_b + nZ(1,0)) - gsl_sf_lnbeta(prior_l2_a, prior_l2_b);
            
            //root
            priorZ += prior_z.at(configZ[i][root]);
            
            // P(X|Z, r*)
            vector<bool> visited = vector<bool>(N,false);
            Z = configZ[i];
            for(int s=S; s<N; s++)
                for(int g=0; g<GG; g++)
                    lambda[g][s].fill(0);
            sum_full_loglik[i] = 0;
            
            for(int g=0; g<GG; g++){
                Update_Tg(g, visited, bpp, false);  //get lambda
                sum_full_loglik[i] += BPP::log_exp_sum(lambda[g][root]);
            }
            
            sum_full_loglik[i] += priorZ;
            
            
        }
        
        VB = BPP::log_exp_sum(conv_to<mat>::from(sum_full_loglik)); //- log(numConf);
        
        VB += log(gsl_ran_gamma_pdf(rc,bpp.cprior_a,bpp.cprior_b));
        if(resZ!=0) VB += log(gsl_ran_gamma_pdf(rn,bpp.nprior_a,bpp.nprior_b));
        
        
        if(verbose) cout << CC <<  " P(X, r): " << VB <<" " << endl;
        
        
        if(resZ==0)  // no or few Z==2, var_nrate[i] < 1e-6
        {
            varc -= pow(meanc,2)/num_mcmc;
            logdet = log(varc) - log(num_mcmc);
            VB += 0.5*logdet + 0.5*log(2*M_PI); //*num_mcmc;
        }else{
            varn -= pow(meann,2)/num_mcmc;
            varc -= pow(meanc,2)/num_mcmc;
            cov -= meann * meanc/num_mcmc ;
            logdet = log(varc * varn - cov * cov) - 2*log(num_mcmc);
            VB += 0.5*logdet + log(2*M_PI); //*num_mcmc;
        }
        
        if(verbose) cout << CC <<  " logdet: " << logdet <<" " << endl;
    }else{
        //chib method
        if(resZ!=0) // full
        {
            VB = VB * ((double)num_mcmc/effective_sum) + log(num_mcmc) + log(2*M_PI) + 1; //effective_sum
            
        }
        else
        {
            VB = VB* ((double)num_mcmc/effective_sum) + log(num_mcmc) + 0.5*(log(2*M_PI) + 1);
            
        }
    }
    
    //chib method
    if(resZ==1) // full
    {
        bpp.log_liks_sgl[CC] = VB;
        
        if(verbose) cout << "P(X,r*|Z*): " << MaxLoglik << " P(X|null): " <<  bpp.log_liks_null[CC] <<" P(X|acc): " << bpp.log_liks_resZ[CC] << " P(X|full): " << bpp.log_liks_sgl[CC] << endl << endl;
        //cout << "P(X,Z*,r*): " << MaxLoglik << " P(X): " << bpp.log_liks_sgl[CC] << endl;

        if(::isnan(bpp.log_liks_sgl[CC])) failure = 1;
    }
    else if(resZ==0)
    {
        bpp.log_liks_null[CC] = VB;
        if(::isnan(bpp.log_liks_null[CC])) failure = 1;

    }else{
        bpp.log_liks_resZ[CC] = VB ;
        if(::isnan(bpp.log_liks_resZ[CC])) failure = 1;
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
            if(Tg[g][s] >= num_base) continue;
            
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
            if(missing[s] || Tg[g][s]>=num_base) continue;
            
            int* p = children2[s];
            lambda[g][s].fill(0);
            
            for(int cc=0;cc<2;cc++)
            {
                int chi = p[cc];
                assert(chi != -1);
                if(missing[chi] || Tg[g][chi]>=num_base) continue;
                
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
//    MH_ratio_g = 0;
//    MH_ratio_l = 0;
//
//
//    for(vector<int>::iterator it = nodes.begin(); it <nodes.end() - 1; it++)
//    {
//        int s = *it;
//        if(Z[parent2[s]] == 0)
//        {
//            //if(log_TM_Int[s](1,0) == 0 || log_TM_Int[s](0,0) == 0) continue;
//            //MH_ratio_g +=  log_TM_Int[s](1,0) * Z[s] + log_TM_Int[s](0,0) * (1 - Z[s]);
//            if(Z[s] == 0) {
//                MH_ratio_g += log_TM_Int[s](0,0);
//                if(fixZ[s] != 1) MH_ratio_l += log_TM_Int[s](0,0);
//            }
//            else if(Z[s] == 1) MH_ratio_g += log_TM_Int[s](1,0);
//            else MH_ratio_l += log_TM_Int[s](2,0);
//
//        }else if(Z[parent2[s]] == 1){
//            if(log_TM_Int[s](2,1) == 0 || log_TM_Int[s](1,1) == 0) continue;
//            MH_ratio_l +=  log_TM_Int[s](2,1) * (Z[s] - 1) + log_TM_Int[s](1,1) * (2 - Z[s]);
//        }
//
//
//    }
    
    // sample g and l rates together
        MH_ratio_g = 0;
    
    
        for(vector<int>::iterator it = nodes.begin(); it <nodes.end() - 1; it++)
        {
            int s = *it;
            if(Z[parent2[s]] < 2)
            {
                MH_ratio_g += log_TM_Int[s](Z[s],Z[parent2[s]]);
            }
            
        }
    
    
    
}


//double BPP_C::prior_Z_subtree(vector<int> & tmpZ)  //whether cal prob_back again, default cal=0
//{
//    double result =0;
//    //vector<vec> log_prob_back_temp = vector<vec> (N, zeros<vec>(3));
//
//    // restore prior of Z (no X, no missing only constraint!)
//    for(vector<int>::iterator it = nodes.begin(); it < nodes.end();it++) // orginally, -1 ...
//    {
//        int s = *it;
//        //if(missing[s] && s<S) continue; //don't update log_prob_back for missing nodes
//        log_prob_back[s].fill(0);
//        if(fixZ[s] == 1)
//        {
//            log_prob_back[s][2] = -INFINITY; //log_prob_back[s][0] = -INFINITY;
//        }
//        //        else if(fixZ[s] ==2)  // no use!
//        //        {
//        //            log_prob_back_temp[s][0] = -INFINITY;  //log_prob_back[s][1] = -INFINITY;
//        //        }
//    }
//
//
//    //message passing from bottom to top
//    //mat tmp = zeros<mat> (3,2);
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
//    log_prob_back[root] += prior_z;
//
//
//
//    //update Z from top to bottom Z[N] = 0 or 1
//    result =log_prob_back[root][tmpZ[root]] - BPP::log_exp_sum(log_prob_back[root]);
//
//    vec log_trans_p, trans_p;
//    for(vector<int>::iterator it =nodes.end()-2;it>=nodes.begin();it--)  //
//    {
//
//        int s = *it;
//        int p = parent2[s];
//
//        if(tmpZ[p]!=2)
//        {
//
//            log_trans_p = log_TM_Int[s].unsafe_col(tmpZ[p]) + log_prob_back[s];
//            result += log_trans_p[tmpZ[s]] - BPP::log_exp_sum(log_trans_p);
//
//            //trans_p = BPP::log_sample(log_trans_p);
//
//
//        }else if(tmpZ[s]!=2){
//            cout <<"prior_Z_subtree error, Z[s]!=2, p: " << p <<"s: "<<s <<endl;
//        }
//
//    }
//
//    return(result);
//    // cout << endl;
//}


//// no use
//vector<double> BPP_C::prior_Z_subtree(vector< vector<int> > & configZ, vector< int > numConfigZ)  // integrate transition priob
//{
//
//    vector<double> result = vector<double>(configZ.size(), 0.0);
//    vector< vector<int> > record_nz_g = vector< vector<int> >(configZ.size(), vector<int> (2, 0));
//    vector< vector<int> > record_nz_l = vector< vector<int> >(configZ.size(), vector<int> (2, 0));
//    vector< vector<int> > record_nz_l2 = vector< vector<int> >(configZ.size(), vector<int> (2, 0));
//
//    // record number of Z transitions for each configuration
//    for(size_t i =0; i < configZ.size(); i++)
//    {
//        mat nZ = zeros(3,2);
//        for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
//        {
//            int p = parent2[*it];
//            if(configZ[i][p] < 2)
//            {
//                nZ(configZ[i][*it],configZ[i][p]) += 1;
//            }
//
//
//        }
//
//        record_nz_g[i][0] = nZ(1,0) + nZ(2,0);
//        record_nz_g[i][1] = nZ(0,0);
//
//
//        //for(vector<int>:: iterator it = upper_c.begin(); it < upper_c.end() -1;it++)
//        for(vector<int>::iterator it = nodes.begin(); it < nodes.end() - 1; it++)
//        {
//            if(fixZ[*it] == 1)
//            {
//                int p = parent2[*it];
//                assert(p!=N);
//                if(configZ[i][p] < 2)
//                {
//                    nZ(configZ[i][*it],configZ[i][p]) -= 1;
//                }
//            }
//        }
//
//
//        record_nz_l[i][0] = nZ(2,1);
//        record_nz_l[i][1] = nZ(1,1);
//
//        record_nz_l2[i][0] = nZ(2,0);
//        record_nz_l2[i][1] = nZ(1,0);
//
//
//
//    }
//
//
//
//    for(size_t i =0; i < configZ.size(); i++)
//    {
//        // get posterior mean of transition rates
//        double grate = (prior_g_a +  record_nz_g[i][0])/(prior_g_a +  record_nz_g[i][0] + prior_g_b +  record_nz_g[i][1]);
//        double lrate = (prior_l_a +  record_nz_l[i][0])/(prior_l_a + prior_l_b +  record_nz_l[i][0] +  record_nz_l[i][1]);
//        double lrate2 = (prior_l2_a +  record_nz_l2[i][0])/(prior_l2_a +  record_nz_l2[i][0] + prior_l2_b +  record_nz_l2[i][1]);
//
//
//        double prior = 0;  // compute prior of transition rates: P(l, g |C) sum over all config of Z
//        for(size_t j =0; j < configZ.size(); j++)
//        {
//             double ll = log((double)numConfigZ[j]/ num_mcmc) + log(gsl_ran_beta_pdf(grate, prior_g_a + record_nz_g[j][0], prior_g_b + record_nz_g[j][1])) + log(gsl_ran_beta_pdf(lrate, prior_l_a +record_nz_l[j][0], prior_l_b + record_nz_l[j][1])) + log(gsl_ran_beta_pdf(lrate2, prior_l2_a + record_nz_l2[j][0], prior_l2_b + record_nz_l2[j][1]));
//            prior += exp(ll);
//        }
//
//        prior = log(prior);
//
//
//        // P(Z | l, g, C)
//        for(vector<int>::iterator it = nodes.begin(); it <nodes.end(); it++)
//        {
//            int s = *it;
//
//            if(fixZ[s] == 1)
//            {
//                log_TM_Int[*it](1,1) = 0;
//                log_TM_Int[*it](2,1) = log(0);
//
//                log_TM_Int[*it](0,0) = log(1 - grate);
//                log_TM_Int[*it](1,0) = log(grate);
//                log_TM_Int[*it](2,0) = log(0);
//
//            }else{
//                log_TM_Int[s](0,0) = log(1 - grate);
//                log_TM_Int[s](1,0) = log(grate) + log(1 - lrate2);
//                log_TM_Int[s](2,0) = log(grate) + log(lrate2);
//
//                double y = 1 - lrate;
//                log_TM_Int[s](1,1) = log(y);
//                log_TM_Int[s](2,1) = log(1-y);
//            }
//
//        }
//
//
//        // compute P(Z|C, rates)
//
//        // restore prior of Z (no X, no missing only constraint!)
//        for(vector<int>::iterator it = nodes.begin(); it < nodes.end();it++) // orginally, -1 ...
//        {
//            int s = *it;
//            //if(missing[s] && s<S) continue; //don't update log_prob_back for missing nodes
//            log_prob_back[s].fill(0);
//            if(fixZ[s] == 1)
//            {
//                log_prob_back[s][2] = -INFINITY; //log_prob_back[s][0] = -INFINITY;
//            }
//
//        }
//
//        //message passing from bottom to top
//        for(vector<int>::iterator it = nodes.begin(); it < nodes.end()-1;it++)
//        {
//            int s = *it;
//            int p = parent2[s];
//
//            log_prob_back[p] += BPP::log_multi2(log_TM_Int[s], log_prob_back[s]);  //matrix * vector
//            //cout << log_prob_back[p] <<endl;
//        }
//
//        log_prob_back[root] += prior_z;
//
//
//
//        //update Z from top to bottom Z[N] = 0 or 1
//        double loglik = log_prob_back[root][configZ[i][root]] - BPP::log_exp_sum(log_prob_back[root]);
//
//        vec log_trans_p, trans_p;
//        for(vector<int>::iterator it =nodes.end()-2;it>=nodes.begin();it--)  //
//        {
//
//            int s = *it;
//            int p = parent2[s];
//
//            if(configZ[i][p]!=2)
//            {
//                log_trans_p = log_TM_Int[s].unsafe_col(configZ[i][p]) + log_prob_back[s];
//                loglik += log_trans_p[configZ[i][s]] - BPP::log_exp_sum(log_trans_p);
//
//            }else if(configZ[i][s]!=2){
//                cout <<"prior_Z_subtree error, Z[s]!=2, p: " << p <<"s: "<<s <<endl;
//            }
//
//        }
//
//       result[i] = loglik + prior - ( log(gsl_ran_beta_pdf(grate, prior_g_a + record_nz_g[i][0], prior_g_b + record_nz_g[i][1])) + log(gsl_ran_beta_pdf(lrate, prior_l_a +record_nz_l[i][0], prior_l_b + record_nz_l[i][1])) + log(gsl_ran_beta_pdf(lrate2, prior_l2_a + record_nz_l2[i][0], prior_l2_b + record_nz_l2[i][1])));
//
//    }
//
//
//    return(result);
//}





void BPP_C::Output_sampling(int iter, string output_path2, BPP &bpp, int resZ){
    
    string outpath_lik = output_path2 +"_mcmc_trace_M" + to_string(resZ) + "_" + to_string(CC) +".txt";
    ofstream out_lik;
    out_lik.precision(8);
    
    #pragma omp critical
    {
      if(iter ==0)
      {
        out_lik.open(outpath_lik.c_str());
        out_lik << "loglik\trate_n\trate_c\tgrate\tlrate\tlrate2\t";
        for(int s =0 ;s<N;s++){  // header: species name
            out_lik<<bpp.nodes_names[s] << "\t";
        }
        out_lik <<endl;
      }else{
        out_lik.open(outpath_lik.c_str(), ios::app);
      }
        
      for(std::size_t i=0; i< trace_loglik.size();i++)
        {
            out_lik<<trace_loglik[i]<<"\t"<<trace_n_rate[i] << "\t"<<trace_c_rate[i]<<"\t"<<trace_g_rate[i]<<"\t"<<trace_l_rate[i]<<"\t"<<trace_l2_rate[i]<<"\t";
            for(int s=0; s<N;s++)
                out_lik<<trace_Z[i][s]<<"\t";
            out_lik <<endl;
        }
    
    out_lik.close();
    }
}


void BPP_C::Output_init(string output_path,string output_path2, BPP &bpp, ofstream & out_Z, int resZ) //resZ no use
{
    
    std::sort(trace_n_rate.begin(), trace_n_rate.end());
    double n_rate =  trace_n_rate[trace_n_rate.size()/2];
    
    std::sort(trace_c_rate.begin(), trace_c_rate.end());
    double c_rate =  trace_c_rate[trace_c_rate.size()/2];
    
    std::sort(trace_g_rate.begin(), trace_g_rate.end());
    double g_rate =  trace_g_rate[trace_g_rate.size()/2];
    
    std::sort(trace_l_rate.begin(), trace_l_rate.end());
    double l_rate =  trace_l_rate[trace_l_rate.size()/2];
    
    std::sort(trace_l2_rate.begin(), trace_l2_rate.end());
    double l2_rate =  trace_l2_rate[trace_l2_rate.size()/2];
    
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
        out_Z << CC << "\t" << n_rate << "\t" << c_rate << "\t" << g_rate << "\t" << l_rate << "\t" << l2_rate;
        for(int s=0; s<N;s++){
            //out_Z <<"\t"<<countZ[s][0];
            for(int k=0;k<4;k++)
                out_Z <<"\t"<<(double)countZ[s][k]/num_mcmc; //?-1
        }
        
        out_Z <<endl;
    }
}









