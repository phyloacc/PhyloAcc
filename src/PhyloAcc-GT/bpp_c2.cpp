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

struct Cmp2
{
    bool operator()(const pair<int, string> &a, const pair<int, string> &b)
    {
        if (a.first == b.first)
        {
            return a.second <= b.second;
        }
        return a.first > b.first;
    }
};

void BPP_C::Output_sampling(int iter, string output_path2, BPP &bpp, int resZ)
{

    string outpath_lik = output_path2 + "_mcmc_trace_M" + to_string(resZ) + "_" + to_string(CC) + ".txt";
    ofstream out_lik;
    out_lik.precision(8);

#pragma omp critical
    {
        if (iter == 0)
        {
            out_lik.open(outpath_lik.c_str());
            out_lik << "iter\tloglik\tindicator\trate_n\trate_c\tpi_A\tGTtop\t";
            for (int s = 0; s < N; s++)
            { // header: species name
                out_lik << bpp.nodes_names[s] << "\t";
            }
            out_lik<<"grate\tlrate";
            out_lik << endl;
        }
        else
        {
            out_lik.open(outpath_lik.c_str(), ios::app);
        }

        for (std::size_t i = 0; i < num_mcmc + num_burn; i++)
        {
            out_lik << iter << "\t" << trace_loglik[i] << "\t" << trace_indicator[i] << "\t" << trace_n_rate[i] << "\t" << trace_c_rate[i] << "\t" << trace_pi[i][0] << "\t" <<trace_GTtopChg[i] << "\t";
            for (int s = 0; s < N; s++)
                out_lik << trace_Z[i][s] << "\t";
            out_lik <<trace_g_rate[i]<<"\t"<<trace_l_rate[i] << endl;
        }

        out_lik.close();
    }
}

void BPP_C::Output_tree(int iter, string outpathG, BPP &bpp, int resZ)
{
    string outpath_Gt = outpathG + "_Gtree" + to_string(resZ) + "_" + to_string(CC) + ".txt";
    ofstream outG;
    outG.precision(8);

    // output top K genetree, K = 5
    set<pair<int, string>, Cmp2> topK;
    set<pair<int, string>>::iterator itbegin;
    for (map<string, simpletree>::iterator it = trace_genetree.begin(); it != trace_genetree.end(); it++) //Han: number of times each tree sampled.
    {
        if (topK.size() < 10) //Han: get the top 5 most freq gene trees
        {
            topK.insert(make_pair(it->second.count, it->first));
        }
        else
        {
            //
            if (topK.rbegin()->first <= it->second.count)
            {
                topK.insert(make_pair(it->second.count, it->first));
                itbegin = topK.begin();
                advance(itbegin, 9);
                int temp_count = itbegin->first;
                for (; itbegin != topK.end(); itbegin++)
                {
                    if (itbegin->first < temp_count)
                        break;
                }

                topK.erase(itbegin, topK.end());
            }
        }
    }

    #pragma omp critical
    if(iter ==0){
        outG.open(outpath_Gt.c_str());
        outG<<"ID\titer\tTreeNum\tprop\tG"<<endl;
    }else{
        outG.open(outpath_Gt.c_str(), ios::app);
    }

    int nTree=0;
    for (itbegin = topK.begin(); itbegin != topK.end(); itbegin++)
    {   
        nTree++;
        string gt = itbegin->second;
        
        //get average branch len
        for (vector<double>::iterator it = trace_genetree[gt].distances.begin(); it != trace_genetree[gt].distances.end(); it++)
        {
            *it /= itbegin->first;
        }
        std::stringstream buffer;
        trace_genetree[gt].printTree(trace_genetree[gt].root, bpp, buffer);
        outG<< CC << "\t" <<iter<<"\t"<<nTree<<"\t"<< (double)itbegin->first / num_mcmc << "\t" << buffer.str() << endl;
    }
    outG.close();
}

void BPP_C::Output_init(string output_path, string output_path2, BPP &bpp, ofstream &out_Z, ofstream &out_tree, int mod_GT)
{

    int mid = num_mcmc / 2;
    std::sort(trace_n_rate.begin() + num_burn, trace_n_rate.begin() + num_mcmc + num_burn);
    double n_rate = trace_n_rate[mid + num_burn]; //Han: medium

    double c_rate = 0;
    for (std::size_t i = num_burn; i < num_mcmc + num_burn; i++)
    {
        c_rate += trace_c_rate[i];
    }
    c_rate /= num_mcmc; //Han mean

    std::sort(trace_g_rate.begin() + num_burn, trace_g_rate.begin() + num_mcmc + num_burn);
    double g_rate = trace_g_rate[mid + num_burn]; //Han: medium

    std::sort(trace_l_rate.begin() + num_burn, trace_l_rate.begin() + num_mcmc + num_burn);
    double l_rate = trace_l_rate[mid + num_burn];

    std::sort(trace_l2_rate.begin() + num_burn, trace_l2_rate.begin() + num_mcmc + num_burn);
    double l2_rate = trace_l2_rate[mid + num_burn];

    vector<vector<int>> countZ = vector<vector<int>>(N, vector<int>(4, 0));
    for (int s = 0; s < N; s++)
    {

        for (std::size_t i = num_burn; i < num_mcmc + num_burn; i++)
        {

            countZ[s][trace_Z[i][s] + 1]++;
        }

        if (missing[s])
        {
            countZ[s][0] = num_mcmc; //set missing s = 1, though Z[s] can be 0/1/2; only missing in upper Z[s] = -1
        }
    }

    // output top K genetree, K = 5
    set<pair<int, string>, Cmp2> topK;
    set<pair<int, string>>::iterator itbegin;
    for (map<string, simpletree>::iterator it = trace_genetree.begin(); it != trace_genetree.end(); it++) //Han: number of times each tree sampled.
    {
        if (topK.size() < 5) //Han: get the top 5 most freq gene trees
        {
            topK.insert(make_pair(it->second.count, it->first));
        }
        else
        {
            //
            if (topK.rbegin()->first <= it->second.count)
            {
                topK.insert(make_pair(it->second.count, it->first));
                itbegin = topK.begin();
                advance(itbegin, 4);
                int temp_count = itbegin->first;
                for (; itbegin != topK.end(); itbegin++)
                {
                    if (itbegin->first < temp_count)
                        break;
                }

                topK.erase(itbegin, topK.end());
            }
        }
    }

    #pragma omp critical
    {
        out_Z << CC << "\t" << n_rate << "\t" << c_rate << "\t" << g_rate << "\t" << l_rate << "\t" << l2_rate;
        for (int s = 0; s < N; s++)
        {
            //out_Z <<"\t"<<countZ[s][0];
            for (int k = 0; k < 4; k++)
                out_Z << "\t" << (double)countZ[s][k] / num_mcmc; //?-1
        }

        out_Z << endl;
    }

    for (itbegin = topK.begin(); itbegin != topK.end(); itbegin++)
    {
        string gt = itbegin->second;
        
        //get average branch len
        if(mod_GT==0){
        for (vector<double>::iterator it = trace_genetree[gt].distances.begin(); it != trace_genetree[gt].distances.end(); it++)
        {
            *it /= itbegin->first;
        }
        }
        std::stringstream buffer;
        trace_genetree[gt].printTree(trace_genetree[gt].root, bpp, buffer);
        #pragma omp critical
        {
            out_tree << CC << "\t" << (double)itbegin->first / num_mcmc << "\t" << buffer.str() << endl;
        }
    }
}


void BPP_C::Output_GTsampling(string output_path2, BPP &bpp, int resZ)
{
    string outpath_mctree = output_path2 + "_trace_genetree_M" + to_string(resZ) + "_" + to_string(CC) + ".txt";
    ofstream out_mctree;
    out_mctree.precision(8);

    #pragma omp critical
    {
        out_mctree.open(outpath_mctree.c_str());
        out_mctree<<"No.\tCount\tGTtopology\tgenetree";
        for(int s=S;s<N;s++){
            out_mctree<<"\t"<<"Node_"<<to_string(s);
        }
        out_mctree<<endl;

        int no = 0;
        for (map<string, simpletree>::iterator it = trace_genetree.begin(); it != trace_genetree.end(); it++)
        {
            no = no + 1;
            string tp = it->first;
            for(vector<double>::iterator it2=trace_genetree[tp].distances.begin(); it2 !=trace_genetree[tp].distances.end();it2++){
                *it2 /= trace_genetree[tp].count;
            } //Get average tree   //If called bppc.Output_init before, then topK5 trees are alr average tree.
            std::stringstream buffer;
            trace_genetree[tp].printTree(trace_genetree[tp].root, bpp, buffer);
            out_mctree << no << "\t" << it->second.count << "\t" << tp << "\t" << buffer.str();
           for(int s=S; s<N; s++){
               out_mctree<<"\t"<<trace_genetree[tp].node_names[s];
           }
            out_mctree<<endl;
        }
    }
    out_mctree.close();
}
