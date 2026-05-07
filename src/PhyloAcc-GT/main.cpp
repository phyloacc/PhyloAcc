//
//  main.cpp
//  PhyloAcc
//
//  Created by hzr on 3/8/16.
//  Copyright © 2016 hzr. All rights reserved.
//

/////////////////////////////////////////////////////////////////

#include <armadillo>

#include <cstdlib>
#include <ctime>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../PhyloAcc-common/newick.h"
#include "../PhyloAcc-common/profile.h"
#include "../PhyloAcc-common/run.h"
#include "bpp.hpp"
#include "bpp_c.hpp"
#include "newick2.h"

using namespace std;
using namespace arma;

/////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    time_t start = time(NULL);

    cout << std::fixed << setprecision(4);

    phyloacc::Config config = phyloacc::LoadConfigForProgram(argc, argv, phyloacc::ProgramKind::GT);
    if(! phyloacc::ValidateOutputDirectory(config))
    {
    	return 1;
    }

    PhyloProf profile = phyloacc::LoadProfile(config);
    phyloacc::DisplayRunSummary(profile, config, phyloacc::ProgramKind::GT);

    PhyloTree phytree = phyloacc::LoadSpeciesTree(config); //Han: .subs_rate contains Q

    PhyloTree_theta tree2;
    if(config.tree_coal_unit !=""){
        tree2 = LoadPhyloTree_theta(config.tree_coal_unit);
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
                        if(phytree.thetas[i]>= config.theta_cutoff){
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
        cout << "Error. Please also input a phylogeny with branch lengths in coalescent units." << endl;
        return 1;
    }

    BPP bpp(0, profile, phytree, config.output_path, config.targetspecies, config.outgroup,
            config.conserve_prop, config.conservegroup, config.ratio0, config.ratio1,
            config.ropt, config.cub, config.nlb, config.nprior_a, config.nprior_b,
            config.cprior_a, config.cprior_b, config.seed,
            config.prep_grate, config.prep_lrate, config.prep_lrate2,
            config.prior_grate_a, config.prior_grate_b, config.prior_lrate_a,
            config.prior_lrate_b, config.prior_lrate2_a, config.prior_lrate2_b,
            config.indel, config.indel2, config.missing_thres, config.sample_indel,
            config.prior_dir_par, config.br_sample_cutoff, config.deepcoal_species);

    bpp.InitMCMC(config.num_burn, config.num_mcmc, config.num_thin);

    phyloacc::RunPaths paths = phyloacc::MakeRunPaths(config);
    phyloacc::OutputBundle outputs = phyloacc::OpenOutputBundle(
        phyloacc::ProgramKind::GT, paths, bpp.nodes_names, config);

    const phyloacc::ModelSpec& m0 = phyloacc::GetModelSpec(phyloacc::ModelId::M0);
    const phyloacc::ModelSpec& m1 = phyloacc::GetModelSpec(phyloacc::ModelId::M1);
    const phyloacc::ModelSpec& m2 = phyloacc::GetModelSpec(phyloacc::ModelId::M2);

    double lrate_prop = 0.5, grate_prop = 0.5;

    vector<int> ids = phyloacc::ResolveElementIds(config, bpp.C, phyloacc::ProgramKind::GT);
    cout << ids.size() << " elements to be computed" << endl;

    if(config.sample_hyper)
    {
        for(int iter  =0; iter<config.num_chain; iter++)
        {
            cout << "Running MCMC chain " << iter +1 << " ..." << endl;
            #pragma omp parallel for schedule (guided) num_threads(config.num_thread)
            for(std::size_t i = 0; i < ids.size(); i++ )
            {
                int c = ids[i];
                bool filter = false;
                bool saw_model_failure = false;
                string completed_models = "none";
                string element_name = phyloacc::ElementName(profile, c);

                try{
                    BPP_C bppc(c, profile, bpp, config.gapchar, config.missing_thres,
                               filter, config.verbose, false, config.consToMis, config.block,
                               config.prune, config.revgap, config.min_length, 0.5, iter);
                    if(filter) {
                        if(config.verbose) cerr << "filter: "<< c <<endl;
                        phyloacc::WriteElementStatus(outputs.status, iter + 1, c, element_name,
                                                     "GT-sample-hyper", "filtered",
                                                     completed_models, "filter");
                        continue;
                    }

                    bppc.initMCMC(0,5,bpp,m2.res_z,config.prune);
                    bppc.Gibbs(0, 4, bpp,outputs.RatePostZ(m2.id),paths.output_prefix,
                               paths.output_prefix2,m2.res_z, true, config.sample_hyper,
                               lrate_prop, grate_prop, false);

                    if(bppc.verbose || bppc.failure) bppc.Output_sampling(iter, paths.output_prefix2, bpp, m2.trace_slot);
                    bppc.Output_init(paths.output_prefix,paths.output_prefix2,bpp,outputs.RatePostZ(m2.id), outputs.Tree(m2.id), bppc.verbose);
                    saw_model_failure = saw_model_failure || bppc.failure;
                    completed_models = m2.suffix;
                    phyloacc::WriteElementStatus(outputs.status, iter + 1, c, element_name,
                                                 "GT-sample-hyper",
                                                 saw_model_failure ? "model_failure" : "ok",
                                                 completed_models, "");

                }catch (exception& e){
                    cout << c << " Standard exception: " << e.what() << endl;
                    phyloacc::WriteElementStatus(outputs.status, iter + 1, c, element_name,
                                                 "GT-sample-hyper", "error",
                                                 completed_models, e.what());
                }
            }
            bpp.sample_hyperparam(iter, ids, outputs.hyper);
            bpp.Output_init0(profile,outputs.likelihood, ids);

        }
    }else if(config.simulate)
    {
        for(std::size_t i = 0; i < ids.size(); i++ )
        {
            int c = ids[i];
            bool filter = false;
            bool saw_model_failure = false;
            string completed_models = "none";
            string element_name = phyloacc::ElementName(profile, c);
            try{
                // accelerate in target species
                BPP_C bppc(c, profile, bpp, config.gapchar, config.missing_thres,
                           filter, config.verbose, false, config.consToMis, config.block,
                           config.prune, config.revgap, config.min_length, 0.5, 0);
                bppc.simulate(bpp, profile, config.gapchar,config.prune);
                saw_model_failure = saw_model_failure || bppc.failure;
                completed_models = "simulate";
                phyloacc::WriteElementStatus(outputs.status, 1, c, element_name,
                                             "GT-simulate", saw_model_failure ? "model_failure" : "ok",
                                             completed_models, "");
            }catch (exception& e){
                cout << c << " Standard exception: " << e.what() << endl;
                phyloacc::WriteElementStatus(outputs.status, 1, c, element_name,
                                             "GT-simulate", "error",
                                             completed_models, e.what());
            }
        }

        bpp.Output_simu(profile, paths.output_prefix, ids.size());
    }else{
        #pragma omp parallel for schedule (guided) num_threads(config.num_thread)
        for(std::size_t i = 0; i < ids.size(); i++ )
        {
            int c = ids[i];
            bool filter = false;
            bool saw_model_failure = false;
            string completed_models = "none";
            string element_name = phyloacc::ElementName(profile, c);

            try{
                BPP_C bppc(c, profile, bpp, config.gapchar, config.missing_thres,
                           filter, config.verbose, config.verboseGT, config.consToMis,
                           config.block, config.prune, config.revgap, config.min_length, 0.5, 0);
                cout<<"element "<<to_string(c)<<", number of base pair="<<to_string(bppc.GG)<<endl;
                if(filter) {
                  if(config.verbose) cerr << "filter: "<< c <<endl;
                  phyloacc::WriteElementStatus(outputs.status, 1, c, element_name,
                                               "GT", "filtered", completed_models, "filter");
                  continue;
                }

                int tot = 0;
                if(bppc.idblk_count==0){
                    double nblk=(double)(bppc.GG - 15)/config.block;
                    int nblk2=(bppc.GG - 15)/config.block;
                    if( (nblk- nblk2)< ((double) 1.0/3.0)){
                        tot=nblk2;
                    }else{
                        tot=nblk2+1;
                    }
                }else{
                    double nblk = (double)(bppc.GG - bppc.idblk_count-15)/config.block;
                    int nblk2 = (bppc.GG - bppc.idblk_count-15)/config.block;
                    if((nblk- nblk2)< ((double) 1.0/3.0)){
                        tot = 1+nblk2;
                    }else{
                        tot = 2+nblk2;
                    }
                }

                cout<<"start null model\n";
                for (int iter = 0; iter <= tot; iter++)
                {
                    bppc.initMCMC(iter, tot, bpp, m0.res_z, config.prune, false);
                    bppc.Gibbs(iter, tot, bpp, outputs.RatePostZ(m0.id), paths.output_prefix,
                               paths.output_prefix2, m0.res_z, true, config.sample_hyper,
                               lrate_prop, grate_prop, config.WL);
                    if (bppc.verbose || bppc.failure)
                    {
                        bppc.Output_sampling(iter, paths.output_prefix2, bpp, m0.trace_slot);
                        bppc.Output_tree(iter, paths.output_prefix2, bpp, m0.trace_slot);
                    }
                }
                bppc.Output_init(paths.output_prefix,paths.output_prefix2,bpp,outputs.RatePostZ(m0.id), outputs.Tree(m0.id), bppc.verbose);
                saw_model_failure = saw_model_failure || bppc.failure;
                completed_models = m0.suffix;

                cout<<"start restricted model\n";
                for (int iter = 0; iter <= tot; iter++)
                {
                    bppc.initMCMC(iter, tot, bpp, m1.res_z, config.prune, false);
                    bppc.Gibbs(iter, tot, bpp, outputs.RatePostZ(m1.id), paths.output_prefix,
                               paths.output_prefix2, m1.res_z, true, config.sample_hyper,
                               lrate_prop, grate_prop, config.WL);
                    if (bppc.verbose || bppc.failure)
                    {
                        bppc.Output_sampling(iter, paths.output_prefix2, bpp, m1.trace_slot);
                        bppc.Output_tree(iter, paths.output_prefix2, bpp, m1.trace_slot);
                    }
                }
                bppc.Output_init(paths.output_prefix,paths.output_prefix2,bpp,outputs.RatePostZ(m1.id), outputs.Tree(m1.id), bppc.verbose);
                saw_model_failure = saw_model_failure || bppc.failure;
                completed_models = m0.suffix + "," + m1.suffix;

                cout<<"start full model\n";
                for (int iter = 0; iter <= tot; iter++)
                {
                    bppc.initMCMC(iter, tot, bpp, m2.res_z, config.prune, false);
                    bppc.Gibbs(iter, tot, bpp, outputs.RatePostZ(m2.id), paths.output_prefix,
                               paths.output_prefix2, m2.res_z, true, config.sample_hyper,
                               lrate_prop, grate_prop, config.WL);
                    if (bppc.verbose || bppc.failure)
                    {
                        bppc.Output_sampling(iter, paths.output_prefix2, bpp, m2.trace_slot);
                        bppc.Output_tree(iter, paths.output_prefix2, bpp, m2.trace_slot);
                    }
                }
                bppc.Output_init(paths.output_prefix,paths.output_prefix2,bpp,outputs.RatePostZ(m2.id), outputs.Tree(m2.id), bppc.verbose);
                saw_model_failure = saw_model_failure || bppc.failure;
                completed_models = m0.suffix + "," + m1.suffix + "," + m2.suffix;

                cout << c << "\t" << bpp.log_liks_WL[0][c] <<"\t" <<  bpp.log_liks_WL[2][c] <<"\t" <<  bpp.log_liks_WL[1][c] <<endl;
                cout<<"\t" << bpp.log_liks_Z[0][c] << "\t" << bpp.log_liks_Z[2][c]<<"\t" << bpp.log_liks_Z[1][c] << endl;
                phyloacc::WriteElementStatus(outputs.status, 1, c, element_name,
                                             "GT", saw_model_failure ? "model_failure" : "ok",
                                             completed_models, "");

            }catch (exception& e){
              cout << c << " Standard exception: " << e.what() << endl;
              phyloacc::WriteElementStatus(outputs.status, 1, c, element_name,
                                           "GT", "error", completed_models, e.what());
            }
      }

      bpp.Output_init(profile,paths.output_prefix, ids);
    }

    outputs.Close();

    cout << endl << endl << "time used:  " << (time(NULL)-start)/60 << " min." << endl << endl;
    return 0;
}
