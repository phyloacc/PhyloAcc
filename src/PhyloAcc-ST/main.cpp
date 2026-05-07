//
//  main.cpp
//  PhyloAcc
//
//  Created by hzr on 3/8/16.
//  Copyright © 2016 hzr. All rights reserved.
//

/////////////////////////////////////////////////////////////////

#include <armadillo>

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

using namespace std;
using namespace arma;

/////////////////////////////////////////////////////////////////
// Main

int main(int argc, char* argv[])
{
    time_t start = time(NULL);

    cout << std::fixed << setprecision(4);

    phyloacc::Config config = phyloacc::LoadConfigForProgram(argc, argv, phyloacc::ProgramKind::ST);
    if (!phyloacc::ValidateOutputDirectory(config))
    {
        return 1;
    }

    PhyloProf profile = phyloacc::LoadProfile(config);
    PhyloTree phytree = phyloacc::LoadSpeciesTree(config);
    phyloacc::DisplayRunSummary(profile, config, phyloacc::ProgramKind::ST);

    BPP bpp(0, profile, phytree, config.output_path, config.targetspecies, config.outgroup,
            config.conserve_prop, config.conservegroup, config.ratio0, config.ratio1,
            config.ropt, config.cub, config.nlb, config.nprior_a, config.nprior_b,
            config.cprior_a, config.cprior_b, config.seed, config.prep_grate,
            config.prep_lrate, config.prep_lrate2, config.prior_grate_a,
            config.prior_grate_b, config.prior_lrate_a, config.prior_lrate_b,
            config.prior_lrate2_a, config.prior_lrate2_b, config.indel,
            config.indel2, config.missing_thres, config.sample_indel);

    bpp.InitMCMC(config.num_burn, config.num_mcmc, config.num_thin);

    phyloacc::RunPaths paths = phyloacc::MakeRunPaths(config);
    phyloacc::OutputBundle outputs = phyloacc::OpenOutputBundle(
        phyloacc::ProgramKind::ST, paths, bpp.nodes_names, config);

    const phyloacc::ModelSpec& m0 = phyloacc::GetModelSpec(phyloacc::ModelId::M0);
    const phyloacc::ModelSpec& m1 = phyloacc::GetModelSpec(phyloacc::ModelId::M1);
    const phyloacc::ModelSpec& m2 = phyloacc::GetModelSpec(phyloacc::ModelId::M2);

    double lrate_prop = 0.5, grate_prop = 0.5;

    vector<int> ids = phyloacc::ResolveElementIds(config, bpp.C, phyloacc::ProgramKind::ST);
    cout << ids.size() << " of elements to be computed" << endl;

    for(int iter =0; iter<config.num_chain; iter++)
    {
        cout << "Running MCMC chain " << iter +1 << " ..." << endl;
        vector<string> rate_rows_m0(ids.size());
        vector<string> rate_rows_m1(ids.size());
        vector<string> rate_rows_m2(ids.size());
        vector<string> status_rows(ids.size());

        #pragma omp parallel for schedule (guided) num_threads(config.num_thread)
        for(std::size_t i = 0; i < ids.size(); i++ )
        {
            int c = ids[i];
            bool filter = false;
            bool saw_model_failure = false;
            string completed_models = "none";
            string element_name = phyloacc::ElementName(profile, c);
            try {
                BPP_C bppc(c, profile, bpp, config.gapchar, config.missing_thres,
                           filter, config.verbose, config.consToMis, config.prune,
                           config.revgap, config.min_length, 1, iter);

                if(filter) {
                    if(config.verbose) cerr << "filter: "<< c <<endl;
                    status_rows[i] = phyloacc::FormatElementStatus(iter + 1, c, element_name,
                                                                   "ST", "filtered",
                                                                   completed_models, "filter");
                    continue;
                }

                if(!config.sample_hyper)
                {
                    bppc.initMCMC(iter,bpp,m0.res_z);
                    bppc.Gibbs(iter,bpp,outputs.RatePostZ(m0.id),paths.output_prefix,
                               paths.output_prefix2,m0.res_z,true,config.sample_hyper,
                               lrate_prop, grate_prop);
                    bppc.Eval2(bpp,m0.res_z);
                    if(bppc.verbose || bppc.failure) bppc.Output_sampling(iter, paths.output_prefix2, bpp, m0.trace_slot);
                    rate_rows_m0[i] = bppc.Output_init_row(bpp, m0.res_z);
                    saw_model_failure = saw_model_failure || bppc.failure;
                    completed_models = m0.suffix;

                    bppc.initMCMC(iter,bpp,m1.res_z);
                    bppc.Gibbs(iter,bpp,outputs.RatePostZ(m1.id),paths.output_prefix,
                               paths.output_prefix2,m1.res_z,true,config.sample_hyper,
                               lrate_prop, grate_prop);
                    bppc.Eval2(bpp,m1.res_z);
                    if(bppc.verbose || bppc.failure) bppc.Output_sampling(iter, paths.output_prefix2, bpp, m1.trace_slot);
                    rate_rows_m1[i] = bppc.Output_init_row(bpp, m1.res_z);
                    saw_model_failure = saw_model_failure || bppc.failure;
                    completed_models = m0.suffix + "," + m1.suffix;
                }

                bppc.initMCMC(iter,bpp,m2.res_z);
                bppc.Gibbs(iter, bpp,outputs.RatePostZ(m2.id),paths.output_prefix,
                           paths.output_prefix2,m2.res_z, true, config.sample_hyper,
                           lrate_prop, grate_prop);
                bppc.Eval2(bpp,m2.res_z);
                if(bppc.verbose || bppc.failure) bppc.Output_sampling(iter, paths.output_prefix2, bpp, m2.trace_slot);
                rate_rows_m2[i] = bppc.Output_init_row(bpp, m2.res_z);
                saw_model_failure = saw_model_failure || bppc.failure;
                completed_models = config.sample_hyper ? m2.suffix : m0.suffix + "," + m1.suffix + "," + m2.suffix;
                status_rows[i] = phyloacc::FormatElementStatus(iter + 1, c, element_name,
                                                               "ST", saw_model_failure ? "model_failure" : "ok",
                                                               completed_models, "");

            }catch (exception& e){
              cout << c << " 1 Standard exception: " << e.what() << endl;
              status_rows[i] = phyloacc::FormatElementStatus(iter + 1, c, element_name,
                                                             "ST", "error", completed_models,
                                                             e.what());
            }
        }

        for(std::size_t i = 0; i < ids.size(); i++)
        {
            outputs.RatePostZ(m0.id) << rate_rows_m0[i];
            outputs.RatePostZ(m1.id) << rate_rows_m1[i];
            outputs.RatePostZ(m2.id) << rate_rows_m2[i];
            outputs.status << status_rows[i];
        }

        try{
          if(config.sample_hyper) {
            bpp.sample_hyperparam(iter, ids, outputs.hyper);
            bpp.Output_init0(profile,outputs.likelihood, ids);
          }else{
            bpp.Output_init(profile,paths.output_prefix, ids);
          }
        }catch (exception& e){
              cout << " 2 Standard exception: " << e.what() << endl;
        }
    }

    outputs.Close();

    cout << endl << endl << "time used:  " << (time(NULL)-start)/60 << " min." << endl << endl;

    return 0;
}

/////////////////////////////////////////////////////////////////
