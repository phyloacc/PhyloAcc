//
//  genetree.cpp
//  PhyloAcc_init3-1
//
//  Created by hzr on 2019/5/19.
//  Copyright © 2019 hzr. All rights reserved.
//
#include <iomanip>
#include <stdio.h>
#include "bpp_c.hpp"
#include "genetree.hpp"

struct Cmp
{
    bool operator ()(const pair<int, double> &a, const pair<int, double> &b)
    {
        
//        if (a.first == b.first)
//      {
//            return a.second < b.second;
//        }

        return a.second <= b.second;
    }
};


void GTree:: copyto(int len, GTree & gtree){
    //gtree.N = N;
    //gtree.S = S;
    gtree.RNG = RNG;
    gtree.GG = len;
    gtree.root = root;
    
    // gene tree
//    gtree.children_gene    = new int[N][2];
//    gtree.parent_gene      = new int[N];
//    gtree.heights_gene   = new double[N];
//    gtree.childID_gene = new int[N];
//    gtree.missing_gene = new bool[N];
    
//    gtree.parent_gene2 = vector<map<int, int>>(N);
//    gtree.temp_coal = vector<vector<int>>(N);
    
    gtree.gene_nodes = gene_nodes;
    gtree.var_br_node = var_br_node;
    
//    gtree.lambda = vector< map<int, vector<mat>> >(N, map<int, vector<mat>>());
//    gtree.Tg = vector < map<int, vector< int> >> (N, map<int, vector<int>>());
    
    for(int s = 0; s< N; s++)
    {
        gtree.children_gene[s][0] = children_gene[s][0];
        gtree.children_gene[s][1] = children_gene[s][1];
        gtree.parent_gene[s] = parent_gene[s];
        gtree.heights_gene[s] = heights_gene[s];
        gtree.childID_gene[s] = childID_gene[s];
        gtree.missing_gene[s] = missing_gene[s];
        gtree.parent_gene2[s] = parent_gene2[s];
        gtree.temp_coal[s] = temp_coal[s];
        
        for(map<int, vector<mat>>::iterator it = lambda[s].begin(); it!=lambda[s].end(); it++)
        {
            gtree.lambda[s][it->first] = vector<mat>(it->second.begin(), it->second.begin()+len);
        }
        
        for(map<int, vector<int>>::iterator it = Tg[s].begin(); it!=Tg[s].end(); it++)
        {
            gtree.Tg[s][it->first] = vector<int>(it->second.begin(), it->second.begin()+len);
           
        }
        
    }
}


void GTree:: copyfrom(int len, GTree & gtree, int numbase){
  
    root = gtree.root;
    gene_nodes = gtree.gene_nodes;
    var_br_node = gtree.var_br_node;
    //cout<<"root and gene_nodes assigned"<<endl;

    for(int s = 0; s< N; s++)
    {
        children_gene[s][0] = gtree.children_gene[s][0];
        children_gene[s][1] = gtree.children_gene[s][1];
        parent_gene[s] = gtree.parent_gene[s];
        heights_gene[s] = gtree.heights_gene[s];
        childID_gene[s] = gtree.childID_gene[s];
        missing_gene[s] = gtree.missing_gene[s];
        parent_gene2[s] = gtree.parent_gene2[s];
        //temp_children[s] = gtree.temp_children[s];
        temp_coal[s] = gtree.temp_coal[s];

        if(s < S) continue;
        lambda[s].clear();
        Tg[s].clear();
        
        
        for(map<int, vector<mat>>::iterator it = gtree.lambda[s].begin(); it!= gtree.lambda[s].end(); it++)
        {
            lambda[s][it->first] = it->second; //vector<mat>(it->second.begin(), it->second.end());
            lambda[s][it->first].resize(GG, zeros<mat>(numbase,3));           
        }
        
        for(map<int, vector<int>>::iterator it = gtree.Tg[s].begin(); it!= gtree.Tg[s].end(); it++)
        {
            Tg[s][it->first] = it->second; //vector<int>(it->second.begin(), it->second.end());
            Tg[s][it->first].resize(GG, -2);          
        }
        
    }
}



/*
void GTree:: initTree(string tree_str, vector<bool> & missing, set<int> & upper, BPP & bpp) // not implement with missing
{
    
    map<string, int> speciesname = map<string, int>();
    for(int s =0; s<S; s++)
    {   //cout<<"s="<<s<<", sp[s]="<<bpp.species_names[s]<<"\t";
        speciesname[bpp.species_names[s]] = s; //Han*: leaf node ordering is the same as species tree.
    }
    PhyloTree tree = LoadPhyloTree_text(tree_str, speciesname);
    
    gene_nodes.clear(); // only used for size~~
    vector<set<int>> temp_children = vector<set<int>>(N);
    for(int i=0; i<N; i++)
    {
        children_gene[i][0] = -1;
        children_gene[i][1] = -1;
        parent_gene[i] = N;
        heights_gene[i] = bpp.heights[i];
        missing_gene[i] = 0;
        parent_gene2[i].clear();
        temp_coal[i].clear();
        
    }

    for(int i=0; i<N; i++)
    {
        int p = -1;
        for(int j=0; j<N; j++)
        {
            if (tree.dag[i][j])
            {
                p++;
                children_gene[i][p] = j;
                childID_gene[j] = p;
                parent_gene[j] = i;
            }
        }   
    }
    
    for(int i = S; i <N; i++)
    {
        heights_gene[i] = heights_gene[children_gene[i][0]] + tree.distances[children_gene[i][0]];
    }

    int node_label=S;
    
    for(int s = 0; s<S; s++) { //s for gene
        if(missing[s]) {
            missing_gene[s] = true;
            if(upper.find(s) != upper.end()) continue;  //missing and in upper
        }
        temp_children[bpp.parent[s]].insert(s);
        temp_children[s].insert(s);
        parent_gene2[s].insert(make_pair(s, s));      
    }
    
    //1-May missing:
    for(int s=S; s<N; s++)
    {
        if(temp_children[s].size() == 0) continue;

        int p=bpp.parent[s];
        double height2=p<N? bpp.heights[p] : INFINITY;
        for(set<int>:: iterator it = temp_children[s].begin(); it!=temp_children[s].end(); it++){
            int cl1=*it;
            int cl2;
            if(heights_gene[parent_gene[cl1]] < height2){
                for(int c=0; c<2; c++){

                }
                gene_nodes.push_back(cl1);
                parent_gene[cl1]=node_label;
                
                node_label++;
            }
        }

    }

    /////old version. no missing
    for(int s = S; s < N; s++) //s for species
    {
        if(temp_children[s].size() == 0) { // for missing
            continue;
        }

        int p = bpp.parent[s];
        double height2 = p < N? bpp.heights[p] : INFINITY;
               
        set<pair<int, double>, Cmp> tocoal;
        vector<int> toErase;
        for(set<int>::iterator it=temp_children[s].begin(); it!=temp_children[s].end(); it++){
            int cl1=*it;
            int node_label = parent_gene[*it];
            if(node_label == N)
            {
                gene_nodes.push_back(*it); // will include root, height2= infinity
                parent_gene2[s].insert(make_pair(cl1, node_label)); // will include root to N
                toErase.push_back(cl1);
                root = cl1;
                
            }else if(heights_gene[node_label] < height2) { //remove "=". if gnode coincide w. spnode, then coalesced in h2 species.
                if(height2-heights_gene[node_label]<1e-8){//precision issue
                    heights_gene[node_label]=height2;
                }else{
                gene_nodes.push_back(cl1);
                parent_gene2[s].insert(make_pair(cl1, node_label));
                toErase.push_back(cl1);
                tocoal.insert(make_pair(node_label, heights_gene[node_label]));
                temp_children[s].insert(node_label);
                }
            }
        }
        if(toErase.size()>0){
            for(int i=0; i<toErase.size();i++) temp_children[s].erase(toErase[i]);
        }
        
        for(set<pair<int, double>, Cmp>::iterator it = tocoal.begin(); it!=tocoal.end(); it++)
        {
            temp_coal[s].push_back(it->first); 
        }
        
        if(p < N)
        {
            for(set<int>:: iterator it = temp_children[s].begin(); it !=temp_children[s].end(); it++ )
            {
                temp_children[p].insert(*it);
                parent_gene2[s].insert(make_pair(*it, *it));
            }
        }        
    }

    for(int s=S; s<N; s++) //s for gene
    {
        if(missing_gene[children_gene[s][0]] && missing_gene[children_gene[s][1]]) missing_gene[s]=true;
    }
}
*/

void GTree:: initTree(string tree_str, BPP & bpp) // not implement with missing
{ 
    map<string, int> speciesname = map<string, int>();
    for(int s =0; s<S; s++)
    {   //cout<<"s="<<s<<", sp[s]="<<bpp.species_names[s]<<"\t";
        speciesname[bpp.species_names[s]] = s; //Han*: leaf node ordering is the same as species tree.
    }
    PhyloTree tree = LoadPhyloTree_text(tree_str, speciesname);
    
    gene_nodes.clear(); // only used for size~~
    vector<set<int>> temp_children = vector<set<int>>(N);
    for(int i=0; i<N; i++)
    {
        children_gene[i][0] = -1;
        children_gene[i][1] = -1;
        parent_gene[i] = N;
        heights_gene[i] = bpp.heights[i];
        missing_gene[i] = 0;
        parent_gene2[i].clear();
        temp_coal[i].clear();
        
    }

    for(int i=0; i<N; i++)
    {
        int p = -1;
        for(int j=0; j<N; j++)
        {
            if (tree.dag[i][j])
            {
                p++;
                children_gene[i][p] = j;
                childID_gene[j] = p;
                parent_gene[j] = i;
            }
        }   
    }
    
    for(int i = S; i <N; i++)
    {
        heights_gene[i] = heights_gene[children_gene[i][0]] + tree.distances[children_gene[i][0]];
    }
    
    for(int s = 0; s<S; s++) {
        temp_children[bpp.parent[s]].insert(s);
        temp_children[s].insert(s);
        parent_gene2[s].insert(make_pair(s, s));      
    }
    
    for(int s = S; s < N; s++)
    {
        if(temp_children[s].size() == 0) { // for missing
            continue;
        }

        int p = bpp.parent[s];
        double height2 = p < N? bpp.heights[p] : INFINITY;
               
        set<pair<int, double>, Cmp> tocoal;
        vector<int> toErase;
        for(set<int>::iterator it=temp_children[s].begin(); it!=temp_children[s].end(); it++){
            int cl1=*it;
            int node_label = parent_gene[*it];
            if(node_label == N)
            {
                gene_nodes.push_back(*it); // will include root, height2= infinity
                parent_gene2[s].insert(make_pair(cl1, node_label)); // will include root to N
                //temp_children[s].erase(temp_children[s].begin());
                toErase.push_back(cl1);
                root = cl1;
                
            }else if(heights_gene[node_label] < height2) { //remove "=". if gnode coincide w. spnode, then coalesced in h2 species.
                if(height2-heights_gene[node_label]<1e-8){//precision issue
                    heights_gene[node_label]=height2;
                }else{
                gene_nodes.push_back(cl1);
                parent_gene2[s].insert(make_pair(cl1, node_label));
                toErase.push_back(cl1);
                tocoal.insert(make_pair(node_label, heights_gene[node_label]));
                temp_children[s].insert(node_label);
                }
            }
        }
        if(toErase.size()>0){
            for(int i=0; i<toErase.size();i++) temp_children[s].erase(toErase[i]);
        }
        
        for(set<pair<int, double>, Cmp>::iterator it = tocoal.begin(); it!=tocoal.end(); it++)
        {
            temp_coal[s].push_back(it->first); 
        }
        
        if(p < N)
        {
            for(set<int>:: iterator it = temp_children[s].begin(); it !=temp_children[s].end(); it++ )
            {
                temp_children[p].insert(*it);
                parent_gene2[s].insert(make_pair(*it, *it));
            }
        }
        
    }
    //std::stringstream buffer;
    //printTree(root, bpp, buffer);
    //cout <<"GT in initTree: "<< buffer.str() << endl;
    
}

void GTree:: initTree(vector<bool> & missing, set<int> & upper, BPP & bpp){
   
    gene_nodes.clear(); // only used for size~~
    vector<vector<int>> temp_children = vector<vector<int>>(N); //what enters into species 
    for(int i=0; i<N; i++)
    {
        children_gene[i][0] = -1;
        children_gene[i][1] = -1;
        parent_gene[i] = -1;
        heights_gene[i] = bpp.heights[i];
        missing_gene[i] = 0;
        parent_gene2[i].clear();
        temp_coal[i].clear();
    
    }
    
    int node_label = S; //Han: create all the internal nodes on GT. node_label: S-(N-1) based on coalescent order.
    
    for(int s = 0; s<S; s++) {
        if(missing[s]) { //leaf: if eg0.8GG missing, True. Non-leaf: both kids missing, True
            missing_gene[s] = true;
            if(upper.find(s) != upper.end()) continue;  //missing and in upper,i.e., outspecies
        }
        
        temp_children[bpp.parent[s]].push_back(s); 
        temp_children[s].push_back(s);
        parent_gene2[s].insert(make_pair(s, s));
        
    }
    
    // sample coalescent time
    for(int ii = S; ii < N; ii++)
    {
        int s = ii; 
        if(temp_children[s].size() == 0) {
            continue;
        }
        
        int p = bpp.parent[s];
        double height2 = p < N? bpp.heights[p] : INFINITY;
        
        std::shuffle ( temp_children[s].begin(),temp_children[s].end(), bpp.twister); //default_random_engine(seed)
        
        size_t nn = temp_children[s].size();
        double current_height = bpp.heights[s];
        size_t i = nn;
        for(; i >1 ; i--)
        {
            double coal = gsl_ran_exponential(RNG, bpp.thetas[s] / ((i-1) * i)) + current_height;
            if(coal <= height2) { //Han: coalesced in current branch
                current_height = coal; //Han: gene node s's height (<0)
                // remove last two children
                int cl1 = temp_children[s][i-1], cl2 = temp_children[s][i-2];
                gene_nodes.push_back(cl1);
                gene_nodes.push_back(cl2);
                parent_gene2[s].insert(make_pair(cl1, node_label)); //(kid, pa)
                parent_gene2[s].insert(make_pair(cl2, node_label));
              
                temp_coal[s].push_back(node_label); //Han? coalesced nodes in s
                
                parent_gene[cl1] = node_label;
                parent_gene[cl2] = node_label;
                
                childID_gene[cl1] = 0;
                childID_gene[cl2] = 1;
                children_gene[node_label][0] =cl1;
                children_gene[node_label][1] =cl2;
                heights_gene[node_label] = coal;
                               
                unsigned long a = gsl_rng_uniform_int(RNG, i-1);
                temp_children[s].insert(temp_children[s].begin() + a,node_label); //nodes entering & in species s.
                //Han://temp_children: 1. coalesced last 2 child genes. Node_label: as a new s and new gene, replacing the coalesnced one
                //Han: always remove the last 2 kids (after shuffle), so node_label can be inserted at position 0 - (i-1).

                if(missing_gene[cl1] && missing_gene[cl2]) missing_gene[node_label] = true;
                node_label ++;//Han: everytime there is a coalescent event, node_label++. It records the 节点 on the gene tree. In total N 节点. Since Root height is Infinity. So all uncoalesced seq will coalesce before reaching root.
                
            }else{
                break;
            }
             
        }
        
        if(p < N)
        {
            //Han: pointer will end after i positions. So only non-coalesced genes will be added.
            for(vector<int>:: iterator it = temp_children[s].begin(); it <temp_children[s].begin() + i; it++ )
            {
                //Han: temp_children includes grand-children. If grand-children genes didnt coalesce b4 pa(s), then it should include to check if coalesce during s-p branch
                temp_children[p].push_back(*it);
                parent_gene2[s].insert(make_pair(*it, *it)); //s: for species node
            }
        }
        
    }
    
    root = node_label - 1;
    gene_nodes.push_back(root); //Han: it records the seq of genes that get coalesced.
    parent_gene2[N-1].insert(make_pair(root, N));
    parent_gene[root] = N;
}

void GTree::initTree_Sptop(vector<bool> & missing, set<int> & upper, BPP & bpp)
{
//missing: from bpp_c, so info from children2. //bpp_c itself does not have kids, pa, heights etc
    vector<vector<int> > temp_children = vector<vector<int> >(N); //what enters into species 
    gene_nodes.clear(); // only used for size~~
    var_br_node.clear(); //gnodes in species that are too short, can change topology
    prob_var_node.clear();
    for(int i=0; i<N; i++)
    {
        children_gene[i][0] = -1; //bpp.children[i][0];
        children_gene[i][1] = -1; //bpp.children[i][1];
        parent_gene[i] = -1;
        heights_gene[i] = bpp.heights[i];  // only need first S
        missing_gene[i] = missing[i];
        parent_gene2[i].clear();
        temp_coal[i].clear();   
    }
    
    for(int i=0; i<S; i++){  
        if(missing_gene[i] && upper.find(i) != upper.end()){
            continue;
        }else{
            parent_gene2[i].insert(make_pair(i, i));
            temp_children[bpp.parent[i]].push_back(i); //index represents species
            temp_children[i].push_back(i);
            //if(bpp.move_br[i]==1) var_br_node.push_back(i); //should go tgr with gene_node. Possible this node is variable, but it has no sib, so not into gene_node.
        } 
    }
    
    int node_label = S;
    for(int s=S; s<N; s++){ //ii for species
        //int s=ii;
        if(temp_children[s].size() == 0) {
            continue;
        }else{
            int p = bpp.parent[s];
            if(temp_children[s].size() == 2){
                gene_nodes.push_back(temp_children[s][0]);
                gene_nodes.push_back(temp_children[s][1]);
                parent_gene2[s].insert(make_pair(temp_children[s][0], node_label)); //(kid, pa)
                parent_gene2[s].insert(make_pair(temp_children[s][1], node_label));
                temp_coal[s].push_back(node_label); 
                parent_gene[temp_children[s][0]] = node_label;
                parent_gene[temp_children[s][1]] = node_label;
                childID_gene[temp_children[s][0]] = 0;
                childID_gene[temp_children[s][1]] = 1;
                children_gene[node_label][0] = temp_children[s][0];
                children_gene[node_label][1] = temp_children[s][1];
                
                if(p < N){
                    heights_gene[node_label] = bpp.heights[s] + bpp.distances[s] * 0.1;
                }else{
                    heights_gene[node_label] = bpp.heights[s] + bpp.distances[temp_children[s][0]]*0.1;
                }

                if (missing_gene[temp_children[s][0]] && missing_gene[temp_children[s][1]]){
                    missing_gene[node_label] = true;
                }else{
                    missing_gene[node_label] = false;
                } 

                //if((missing_gene[node_label]==false) && (bpp.move_br[s]==1)) var_br_node.push_back(node_label);
                if((missing_gene[node_label]==false) && (bpp.move_br[s]==1) && (p<N)){ //now move_br present short br, not kids of it
                    var_br_node.push_back(temp_children[s][0]); //if br s is short (sp), then its kid's gene can have deep coalesce. 
                    var_br_node.push_back(temp_children[s][1]);
                    double tmp_pr = (double) 1.0/max(0.0005,bpp.heights[p]-bpp.heights[s]);
                    prob_var_node.push_back(tmp_pr);
                    prob_var_node.push_back(tmp_pr);
                }

                temp_children[s].clear();
                temp_children[s].push_back(node_label);
                node_label++;
            }
            if (p < N)
            {
                if(temp_children[s].size() !=1) throw runtime_error("initTree_Sptop error!");
                temp_children[p].push_back(temp_children[s][0]); // only one element in temp_children[s]
                parent_gene2[s].insert(make_pair(temp_children[s][0],temp_children[s][0])); //for gnode with 1kid missing, the other kid coming in and out of this node, and passed to node p
            }
        }
    }
    
    root = node_label - 1;
    gene_nodes.push_back(root); //Han: it records the seq of genes that get coalesced.
    parent_gene2[N-1].insert(make_pair(root, N));
    parent_gene[root] = N;

    //printSptree(bpp, 1);
}

/*** compute the prior probability of the tree ***/
double GTree::priorTree(BPP & bpp)
{
    double logP = 0;
   
    for(int s = S; s<N; s++)
    {
        int nn = parent_gene2[s].size() - temp_coal[s].size(); //Han: num of genes entering sp s
        if(nn == 0) continue;
        double current_height = bpp.heights[s];

        for(vector<int>::iterator it = temp_coal[s].begin(); it != temp_coal[s].end(); it++)
        {
            logP -= (heights_gene[*it] - current_height) * (nn-1) * nn/ bpp.thetas[s];
            logP += log(2) - log(bpp.thetas[s]); //2/theta * exp(-k(k-1)/\theta * t); theta = 4N\mu
            nn--;
            current_height = heights_gene[*it];
        }
        if(nn > 1 && bpp.parent[s]!=N)
        {
            logP -= (bpp.heights[bpp.parent[s]] - current_height) * (nn-1) * nn / bpp.thetas[s] ;
        }
    }
    
    return(logP);
}

/** initial Tg missing **/
void GTree::InitTg(int len, BPP & bpp, vector<int> & nodes, int start)
{
    for(vector<int>::iterator it = nodes.begin();it<nodes.end();it++)  //
    {
        int ss = *it;
        if(parent_gene2[ss].size() ==0 ) continue;  // for missing
        
        vector<int> curr_coal = temp_coal[ss]; // temp_coal needs to be from bottom to top
        
        //for leaf node, below for loop will be skipped.
        //for: for each coalesced node in sp ss.
        for(vector<int>::iterator it2 = curr_coal.begin(); it2 < curr_coal.end(); it2++)
        {
            int gs = *it2;
            
            for(int g = start; g < len; g++)
            {
                int ct = 0; // count Tg missing child
                for(int cc=0;cc<2;cc++)
                {
                    int chi = children_gene[gs][cc];                   
                    if(missing_gene[chi] || Tg[ss].at(chi)[g] >= bpp.num_base ) {
                        ct ++;
                        //continue; //redundant
                    }
                }
                if(ct == 2)  {
                    Tg[ss][gs][g] = bpp.num_base;
                }
                else{
                    Tg[ss][gs][g] = -2;
                }
            }
            
        }
        
        int sp = bpp.parent[ss];
        if(sp < N)  // send message to parent //Han: for not coalesced nodes (passing on to sp branch), its seq bp missing status is inherited from ss when passing to sp.
        {
            for(map<int, int>::iterator it = parent_gene2[ss].begin(); it!= parent_gene2[ss].end(); it++)
            {
                if(it->first != it->second) continue;
                int chi = it->first;
                
                //rewrite for speed
                if(missing_gene[chi])
                {
                    for(int g = start; g < len; g++) Tg[sp].at(chi)[g] = Tg[ss].at(chi)[g];
                }else{
                    
                    for(int g = start; g < len; g++)
                    {
                        if(Tg[ss].at(chi)[g] >= bpp.num_base)
                        {
                            Tg[sp].at(chi)[g] = Tg[ss].at(chi)[g];

                        }else{
                            Tg[sp].at(chi)[g] = -2;

                        }
                    }
                }
            }
        }
    }   
}

void GTree::getGeneNodes(int s, int & minss) // DFS, reorder the children, min species first
{
    
    if(s == -1) {
        minss = -1;
        return;
    }
    int minss1, minss2;
    getGeneNodes(children_gene[s][0], minss1);
    getGeneNodes(children_gene[s][1], minss2);
    if(minss1 == -1)
    {
        minss = s;
    }else if(minss2 < minss1) {
        int temp = children_gene[s][1];
        children_gene[s][1] = children_gene[s][0];
        children_gene[s][0] = temp;
        
        minss = minss2;
    }else{
        minss = minss1;
    }
    
}

// relabel and get branch length by DFS using children_gene (parent_gene for branch length)
void GTree::printTree2(int s, int & currentS, vector<int> & parent_relabel, vector<double> & branchlen)
{
    if (children_gene[s][0] == -1){
        branchlen[s] = heights_gene[parent_gene[s]] - heights_gene[s];
    }else{
        vector<int> currentChild(2);
        for(int i =0;i <2; i++){
            int child = children_gene[s][i];
            printTree2(child, currentS, parent_relabel, branchlen);
            if(child < S){
                currentChild[i] = child;
            }else{
                currentChild[i] = currentS;
                currentS++;
            }          
        }        
        for(int i =0;i <2; i++){
            parent_relabel[currentChild[i]] = currentS;
        }       
        if(parent_gene[s] < N){
             branchlen[currentS] = heights_gene[parent_gene[s]] - heights_gene[s];
        }else{
            branchlen[currentS] = 0;
            parent_relabel[currentS] = N;
        }       
    }
}


void GTree::printTree(int s, BPP& bpp, std::stringstream & buffer)
{
    if (children_gene[s][0] == -1)
    {
        buffer << bpp.species_names[s] << ":"<< heights_gene[parent_gene[s]] - heights_gene[s];
    }
    else
    {
        buffer << "(";
        for(int i =0;i <2; i++)
        {
            int child = children_gene[s][i];
            
            printTree(child, bpp, buffer);
            if(i==0) buffer << ",";
            
        }

        if(parent_gene[s] < N)
        {
            buffer << "):" << heights_gene[parent_gene[s]] - heights_gene[s];
        }else{
            buffer << ");" ;
        }
    }
}

 /** remove the selected branch from gene tree, return the species where the branch coalescent **/
int GTree::Remove_branch(int branch, int ss, BPP & bpp) //branch: branch to be removed on GT; ss: the species this branch is in on Sp tree: i.e. the node this branch leads to on GT is within this ss on sp tree.
{   
    int sp = bpp.parent[ss]; 
    int branchp = parent_gene[branch], branchpp = parent_gene[branchp];
    int branchsib = children_gene[branchp][1 - childID_gene[branch]];
     
    parent_gene[branchsib] = branchpp; //branch is removed from GT, so sib directly links to grandpa.
    if(branchpp != N)
    {
        children_gene[branchpp][childID_gene[branchp]] = branchsib;
        childID_gene[branchsib] = childID_gene[branchp];
    }else{
        root = branchsib; 
    }
    
    bool findcol = false, findp = false;
    int start_ss = ss;  // start species for updating lambda
    
    if(ss < S) ss = sp; //0517, will not be triggered
    int ss_temp =  ss;
    
    while(ss_temp != N)
    {
           int sp_temp=bpp.parent[ss_temp];
           // remove branch above ss, keep coalescent at branch
           if(!findp && ss_temp!=ss) //Han: first run, False. //Han:
           {
               //Han: given branch didn't coalesce in the previous branch, in the next level, if will surely be evaluted.
                if(parent_gene2[ss_temp].find(branch) != parent_gene2[ss_temp].end())
                {
                   lambda[ss_temp].erase(branch);
                   Tg[ss_temp].erase(branch);
                   parent_gene2[ss_temp].erase(branch);
               }
           }
           // remove branchp from ss and set branchp to branchsib
           //Han: branchp is the key. If branch coalesced in ss_temp's branch.(even if branchp not coalesced, pair(branchp, branchp) is in)
           if(parent_gene2[ss_temp].find(branchp) != parent_gene2[ss_temp].end())
           {
               findp = true;
               
               if(!findcol)
               {
                   //Han: return the 1st occurence of branchp in temp_coal[ss]. If not found, return last position of temp_coal[ss]
                   vector<int>::iterator it = find(temp_coal[ss_temp].begin(), temp_coal[ss_temp].end(), branchp);
                   if(it != temp_coal[ss_temp].end())//it will always find branchp, as outer if statement means branchp will be in temp_coal.
                   {
                       temp_coal[ss_temp].erase(it);//Han: if branch coalesced with its sibling to form branchp, then remove branchp, as branch is removed.
                       findcol = true;
                       start_ss = ss_temp;//Han: if branch=leaf, this will update start_ss, o.w. no change
                   }
               }

               lambda[ss_temp].erase(branchp);
               Tg[ss_temp].erase(branchp);
               
               if(parent_gene2[ss_temp][branchp] == branchp) {//Han: if branchp didn't coalesce in ss_temp branch
                   parent_gene2[ss_temp][branchsib] = branchsib; //Han: since branchp is removed, branchsib remains till end of ss_temp

                   lambda[sp_temp].insert(make_pair(branchsib, vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
                   Tg[sp_temp].insert(make_pair(branchsib, vector<int>(GG, -2)));

               }
               else parent_gene2[ss_temp][branchsib] = parent_gene2[ss_temp][branchp]; //Han: if branchp coalesced, sibling's gene parent becomes branchp's gene parent
               
               parent_gene2[ss_temp].erase(branchp);
           }else if(findp) break;
           ss_temp = sp_temp;
    }
    
    return(start_ss);
}

void GTree::Graft(int branch, int ss, int target, double coal, BPP& bpp, double rate, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv)
{
    int branchp = parent_gene[branch];
    int sibid = 1 - childID_gene[branch];
    int it2p = parent_gene[target]; //parent gene of the sibling that will become branchpp
    vector<mat> log_TMs = vector<mat>(2);
    bool findp2 = false;
    
    parent_gene[branchp] = it2p;
    parent_gene[target] = branchp;
    
    children_gene[branchp][sibid] = target;
    if(it2p < N) {
        int cc = childID_gene[target];
        children_gene[it2p][cc] = branchp;
    }else{
        root = branchp;
        
    }
    
    //cout<<"Graft start: ss="<<ss<<", branchp="<<branchp<<", it2p="<<it2p<<", target="<<target<<endl;
    childID_gene[branchp] = childID_gene[target]; //branchp's childID in branchpp inherited from target
    childID_gene[target] = sibid;
    
    vector<int>::iterator tempit = temp_coal[ss].begin();
    for(; tempit != temp_coal[ss].end(); tempit++)
    {
        if(heights_gene[*tempit] > coal)
        {
            temp_coal[ss].insert(tempit,branchp);// Han: insert branchp. temp_coal: small height front.
            break;
        }
    }
    if(tempit == temp_coal[ss].end())
    {
        temp_coal[ss].push_back(branchp);
    }
    
    if(parent_gene2[ss][target] == it2p) {//Han: if previously target coalesced in ss, now branchp will coalesce into it2p, forming (branchp, it2p).
        findp2 = true;
        parent_gene2[ss][branchp] = it2p;
    }else{
        parent_gene2[ss][branchp] = branchp;//Han: if previously target didn't coalesce in ss, branchp will also not, and pair will be (branchp, branchp)
    }
    
    parent_gene2[ss][target] = branchp;
    parent_gene2[ss][branch] = branchp;
    
    //update lambda[ss][branchp], Tg[ss][branchp],
    lambda[ss].insert(make_pair(branchp, vector<mat>(GG,zeros<mat>(bpp.num_base, 3))));
    Tg[ss].insert(make_pair(branchp, vector<int>(GG,-2)));
    
    for(int cc=0;cc<2;cc++)
    {
        int chi = children_gene[branchp][cc];
        double height_child = heights_gene[chi] > bpp.heights[ss] ? heights_gene[chi]:bpp.heights[ss];//Han: branch_height definitely >height[ss]; but its sibling may not. it may come from branches below ss.
        double lentoCal=heights_gene[branchp] - height_child;
       // if(lentoCal<0){lentoCal=1e-9;}
        log_TMs[cc] = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv);//Han: first argument gives "t" in paper, i.e. branch_length in gene_tree between two gene nodes //Transition prob
    }
    
    //cout<<"ss="<<ss<<endl;
    for(int g = 0; g < GG; g++)
    {   
        int ct = 0; // count Tg missing child
        for(int cc=0;cc<2;cc++)
        {   
            int chi = children_gene[branchp][cc];
            //cout<<"chi="<<chi<<endl;
            if(missing_gene[chi] || Tg[ss].at(chi)[g] >= bpp.num_base ) {
                ct ++;
                continue;
            }
            
            lambda[ss][branchp][g].col(cc) = BPP::log_multi(log_TMs[cc], lambda[ss][chi][g].col(2));//Han: lambda[ss][chi] alr computed. Likelihood from gene node sibling, or gene node branch, all unchanged.
            
        }
        lambda[ss][branchp][g].col(2) = lambda[ss][branchp][g].col(1) + lambda[ss][branchp][g].col(0);
        if(ct == 2){
            Tg[ss][branchp][g] = bpp.num_base;
        } 
        
        
    }
    
    if(missing_gene[branch] && missing_gene[target])
    {
        missing_gene[branchp] = true;
    }else{
        missing_gene[branchp] = false;
    } 
    //int ss_temp = ss;
    // from ss to it2p
    while(!findp2 && ss < N)
    {
        ss = bpp.parent[ss];
        lambda[ss].insert(make_pair(branchp, vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
        Tg[ss].insert(make_pair(branchp, vector<int>(GG,-2)));
        lambda[ss].erase(target);
        Tg[ss].erase(target);
        // update parent_gene2
        if(parent_gene2[ss][target] == target)
        {
            parent_gene2[ss][branchp] = branchp;
            
        }else{
            parent_gene2[ss][branchp] = it2p;
            findp2 = true;
        }
        parent_gene2[ss].erase(target);
        
    }
}

bool GTree::Sample_tree(int indicator, BPP & bpp, vector<int> & Z, double n_rate, double c_rate, double consToMis, double nconsToMis, vector<int> lens, vec & loglik,  vec & logp_Z, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv, vec& log_pi)
{
    
    bool iscoal = false;
    bool topaccept = false;
    // randomly select a branch can't be root and can't
    int branch = gsl_rng_uniform_int(RNG, var_br_node.size() - 1); // gene nodes may not contain all 1-S
    branch = gene_nodes[var_br_node[branch]];
    if(branch == root) branch = gene_nodes[gene_nodes.size() - 1];
    
    // get the lineage on species tree
    int ss = branch;
    while(children_gene[ss][0]!=-1)
    {
        ss = children_gene[ss][0];
    }
    
    int sp = bpp.parent[ss];
    double current_height = heights_gene[branch];
    
    while(sp < N && bpp.heights[sp] < current_height) //coalescent strictly higher than speciation, 0517, should it be '<' ? for leaves
    {
        //if(heights_gene[parent_gene[ss]] > current_height) tocoal.insert(s);
        ss = sp ;
        sp = bpp.parent[ss];
    }
    
    int branchp = parent_gene[branch];
    int branchsib = children_gene[parent_gene[branch]][1 - childID_gene[branch]];
    
    int start_ss = Remove_branch(branch, ss, bpp);
    //Han*
    //Update_Lambda(start_ss, branchsib, bpp, Z, n_rate, c_rate); // also update missing pattern of Tg
    Update_Lambda(start_ss, branchsib, bpp, Z, n_rate, c_rate,c_eigenvec, c_eigenval, c_eigeninv); 
    
    // graft the select branch
    while(ss != N)
    {
        // get children of ss not yet coalescent at gene_nodes[branch]
        set<pair<int, double>, Cmp> tocoal;
        double height2 = sp < N? bpp.heights[sp] : INFINITY;
        int nn = 0;
        
        double rate = 1;
        if(Z[ss] == 2)
        {
            rate = n_rate;
        }else if(Z[ss] == 1)
        {
            rate = c_rate;
        }
        
        for(map<int, int>::iterator it = parent_gene2[ss].begin(); it != parent_gene2[ss].end(); it++)
        {
            int gs = it->first;
            if(gs == branch) continue;  // kept in temp_children (leaves) and temp_coal for ss
            
            if(parent_gene[gs] == N || heights_gene[parent_gene[gs]] > current_height)
            {
                tocoal.insert(make_pair(gs, heights_gene[gs]));
                if(heights_gene[gs] <= current_height) nn++; // ?0517
            }
            
        }
        
        if(tocoal.size() >0 )
        {
            set<pair<int, double>>::iterator it = tocoal.begin();
            advance(it,nn-1);
            double current_height2;
            //Move
            while(it != tocoal.end())
            {
                advance(it,1);
                if(it != tocoal.end()) {
                    current_height2 = it->second;
                }
                else current_height2 = height2;
                
                // sample coalescent time
                double coal = gsl_ran_exponential(RNG, bpp.thetas[ss] / (2*nn)) + current_height; //Han: added 2* factor in.
                if(coal <= current_height2) {
                    //getcoal = true;
                    heights_gene[branchp] = coal;
                    
                    set<pair<int, double>>::const_iterator it2(tocoal.begin());
                    int a = gsl_rng_uniform_int(RNG, nn);
                    advance(it2,a);
                    
                    Graft(branch, ss, it2->first, coal, bpp, rate, c_eigenvec, c_eigenval, c_eigeninv);
                    iscoal = true;
                    break;
                }
                
                if(it == tocoal.end()) break;

                int cl1=children_gene[it->first][0];
                int cl2=children_gene[it->first][1]; 
                set<pair<int,double>>::iterator it_erase=tocoal.begin();              
                while(it_erase !=tocoal.end()){
                    int itelem=it_erase->first;
                    if(itelem==cl1 || itelem==cl2){
                        set<pair<int,double>>::iterator it3=it_erase;
                        advance(it_erase,1);
                        tocoal.erase(it3);
                    }else{
                        advance(it_erase,1);
                    }
                }
                nn--;
                current_height = current_height2;
                
            }
        }
        
        if(iscoal) break;
        
        double height_child = heights_gene[branch] > bpp.heights[ss] ? heights_gene[branch]:bpp.heights[ss];
        mat log_TM = bpp.getlogTMc(bpp.heights[sp] - height_child, rate, c_eigenvec, c_eigenval, c_eigeninv);
        
        lambda[sp].insert(make_pair(branch, vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
        Tg[sp].insert(make_pair(branch, vector<int>(GG,-2)));
        
        for(int g = 0; g < GG; g++)
        {
            if(missing_gene[branch] || Tg[ss][branch][g] >= bpp.num_base) {
                Tg[sp][branch][g] = bpp.num_base;
            }else{
                
                lambda[sp][branch][g].col(2) = BPP::log_multi(log_TM, lambda[ss].at(branch)[g].col(2));
            }
        }
        
        
        parent_gene2[ss][branch] = branch;
        current_height = bpp.heights[sp];
        
        ss = sp;
        sp = bpp.parent[ss];
        
    }
    
    Update_Lambda(ss, branchp, bpp, Z, n_rate, c_rate, c_eigenvec, c_eigenval, c_eigeninv);
    
    // MH move
    double r = gsl_rng_uniform(RNG);
    vector<double> loglik_new = vector<double>(2,0);
    for(int g = 0; g< lens[indicator]; g++)
    {
        double temp = BPP::log_exp_sum(lambda[N-1][root][g].col(2) + log_pi);
        loglik_new[0] += temp;
        if(g < lens[1]) loglik_new[1] += temp;
    }
    
    double logp_Z_new = get_logpZ1(bpp,Z, consToMis, nconsToMis); 
    double delta_loglik = loglik_new[indicator] - loglik[indicator];
    
    //cout << "sample tree: "<< delta_loglik <<", " << logp_Z_new << ", " << logp_Z[indicator] << endl;
    if(logp_Z[indicator]!=0) // for iter =0 , indicator =1
    {
        delta_loglik += logp_Z_new - logp_Z[indicator];
    }

    if(log(r) < delta_loglik) // accept
    {
        topaccept = true;
        loglik = loglik_new;
        if(logp_Z[1]==0) // for iter = 0
        {
            logp_Z[0] = logp_Z_new;
        }else{
            logp_Z[0] = logp_Z[1] = logp_Z_new;
        }
    }   
    return(topaccept);
}

double GTree::get_logpZ(vector<int> & Z, double consToMis, double nconsToMis) //assume gene loss happen at the root side of the branch
{
    double logp_Z =0;
    for(int s = S; s < N; s++) //N-1
    {
        for(vector<int>::iterator it = temp_coal[s].begin(); it != temp_coal[s].end(); it++)
        {
            int ss=*it;
            if(missing_gene[ss]){ //if ss is missing, both its kids must be missing. P(kid[0/1] missing | ss missing)=1
                continue;
            }
            for(int cc =0; cc<2; cc++)
            {
                int cl = children_gene[ss][cc];
                if(missing_gene[cl]) //using cl & s. so loss is at root side.
                {
                    logp_Z += log(consToMis) * (Z[s] == 1) + log(nconsToMis) * (Z[s] != 1);
                }else{
                    logp_Z += log(1 - consToMis) * (Z[s] == 1) + log(1 - nconsToMis) * (Z[s] != 1);
                }
            }
        }
    }
    //root if missing, not possible. this element will be ignored.
    //logp_Z+=log(1-consToMis) * (Z[N-1]==1)+ log(1 - nconsToMis) * (Z[N-1] != 1);
    return(logp_Z);
}

/*
double GTree::get_logpZ1(BPP & bpp, vector<int> & Z, double consToMis, double nconsToMis) //assume gene loss happen at the bottom side of the branch
{
    double logp_Z =0;
    for(int s = 0; s < N; s++) //N-1
    {
        for(map<int,int>::iterator it = parent_gene2[s].begin(); it != parent_gene2[s].end(); it++)
        {
            int ss= it->first;
            int sp = parent_gene[ss];
            if(sp != N &&  heights_gene[ss] < bpp.heights[s] && missing_gene[sp]){ //if sp is missing, both its kids must be missing. P(kid[0/1] missing | ss missing)=1
                continue;
            }
           
            if(missing_gene[ss])
            {
                logp_Z += log(consToMis) * (Z[s] == 1) + log(nconsToMis) * (Z[s] != 1);
            }else{
                logp_Z += log(1 - consToMis) * (Z[s] == 1) + log(1 - nconsToMis) * (Z[s] != 1);
            }
            
        }
    }
    return(logp_Z);
}
*/

double GTree::get_logpZ1(BPP & bpp, vector<int> & Z, double consToMis, double nconsToMis) //assume gene loss happen at bottom side of the branch
{
    double logp_Z =0;

    for(int s = 0; s<S; s++){ //tip node gene index = specie index
        if (missing_gene[s])
        {
            logp_Z += log(consToMis) * (Z[s] == 1) + log(nconsToMis) * (Z[s] != 1);
        }
        else
        {
            logp_Z += log(1 - consToMis) * (Z[s] == 1) + log(1 - nconsToMis) * (Z[s] != 1);
        }
    }

    for(int s = S; s < N; s++) //N-1
    {
        for(map<int,int>::iterator it = parent_gene2[s].begin(); it != parent_gene2[s].end(); it++)
        {
            int ss= it -> first;
            if(heights_gene[ss] >= bpp.heights[s]){
                if (missing_gene[ss])
                {
                    logp_Z += log(consToMis) * (Z[s] == 1) + log(nconsToMis) * (Z[s] != 1);
                }
                else
                {
                    logp_Z += log(1 - consToMis) * (Z[s] == 1) + log(1 - nconsToMis) * (Z[s] != 1);
                }
            }
        }
    }
    return(logp_Z);
}

/*
double GTree::get_logpZ1(BPP & bpp, vector<int> & Z, double consToMis, double nconsToMis) //assume gene loss happen at top side of the branch
{
    double logp_Z =0;

    for(int s = 0; s < N; s++) //N-1
    {
        for(map<int,int>::iterator it = parent_gene2[s].begin(); it != parent_gene2[s].end(); it++)
        {
            int ss= it->first;
            int ss2 = it -> second;
            int sp = parent_gene[ss];
            if(sp !=N && missing_gene[sp]){ //if pa is missing, kids must be missing P(missing(kid) | pa) =1.
                continue;
            }else if(sp !=N && !missing_gene[sp] && missing_gene[ss] && ss==ss2){ //if parent not missing, self missing. Assume missing at Top, ie immediate after coal gnode sp.
                continue;
            }
            
            //nodes/brpart counted towards logP(MissingZ) includes:
            //1. root (cannot missing by construction)
            //2. gene not missing, then P(not Missing | Z) needed for all species along its lineage.
            //3. gene missing, then P(Missing | Z_sp=specie for gnode(pa(gene)))
            if(missing_gene[ss])
            {
                logp_Z += log(consToMis) * (Z[s] == 1) + log(nconsToMis) * (Z[s] != 1);
            }else{
                logp_Z += log(1 - consToMis) * (Z[s] == 1) + log(1 - nconsToMis) * (Z[s] != 1);
            }
            
        }
    }
    return(logp_Z);
}
*/

void GTree:: Update_Lambda(int start_ss, int branchsib, BPP& bpp, vector<int> & Z, double n_rate, double c_rate, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv)  
{
    
    // 1. sending the lambda msg from leaves bottom up through the network
    int ss = start_ss;  //where branchp locates
    int current_s = branchsib; //in Sample_Tree, branchsib can be branchp
    while(ss < N)
    {
        int current_p = parent_gene2[ss][current_s];
        
        double rate = 1;
        if(Z[ss] == 2)
        {
            rate = n_rate;
        }else if(Z[ss] == 1)
        {
            rate = c_rate;
        }
               
        while(current_s != current_p && current_p != N)
        {
            int cc = childID_gene[current_s];  // already modified  childID_gene[branchp] = childID_gene[branchsib]
            int otherc = children_gene[current_p][1 - cc];
            double height_child = heights_gene[current_s] > bpp.heights[ss] ? heights_gene[current_s]:bpp.heights[ss];
            double lentoCal=heights_gene[current_p] - height_child;
            //if(lentoCal<0){lentoCal=1e-9;}
            mat log_TM = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv);
           
            for(int g = 0; g < GG; g++)
            {
                int ct = 0;
                if(missing_gene[current_s] || Tg[ss].at(current_s)[g] >= bpp.num_base)
                {
                    lambda[ss].at(current_p)[g].col(cc).fill(0);
                    ct ++;
                }else{
                    lambda[ss][current_p][g].col(cc) = BPP::log_multi(log_TM, lambda[ss].at(current_s)[g].col(2));
                }
                
                if(missing_gene[otherc] || Tg[ss].at(otherc)[g]  >= bpp.num_base)
                {
                    ct ++ ;
                }
                   
                if(ct == 2)
                {
                    Tg[ss].at(current_p)[g] = bpp.num_base;
                }else{
                    Tg[ss].at(current_p)[g] = -2;
                }
                                       
                lambda[ss][current_p][g].col(2) = lambda[ss][current_p][g].col(1) + lambda[ss][current_p][g].col(0);               
            }
            
            if(missing_gene[children_gene[current_p][0]] && missing_gene[children_gene[current_p][1]])
            {
                missing_gene[current_p] = true;
            }else{
                missing_gene[current_p] = false;
            }            
            current_s = current_p;
            current_p = parent_gene2[ss][current_p];
        }
                
        int sp = bpp.parent[ss];
        if(sp < N)  // send message to parent
        {
            
            double height_child = heights_gene[current_s] > bpp.heights[ss] ? heights_gene[current_s]:bpp.heights[ss];
            double lentoCal=bpp.heights[sp] - height_child;
            if(lentoCal<0){lentoCal=1e-9;}
            mat log_TM = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv);
            
            if(missing_gene[current_s])
            {
                Tg[sp].at(current_p) = Tg[ss].at(current_p);
                for(int g = 0; g < GG; g++)
                {
                    //Tg[sp].at(current_p)[g] = Tg[ss].at(current_p)[g]; //11Feb-21 Debug //cur_s=cur_p
                    lambda[sp][current_p][g].fill(0); //= zeros<mat>(bpp.num_base, 3);
                }
            }else{
                for(int g = 0; g < GG; g++)
                {
                    if(Tg[ss].at(current_s)[g] >= bpp.num_base)
                    {
                        Tg[sp].at(current_p)[g] = Tg[ss].at(current_p)[g]; //11Feb-21 Debug //cur_s=cur_p
                        lambda[sp][current_p][g].fill(0); //= zeros<mat>(bpp.num_base, 3);
                    }else{
                        lambda[sp][current_p][g].col(2) = BPP::log_multi(log_TM, lambda[ss].at(current_s)[g].col(2));
                        Tg[sp].at(current_p)[g] = -2;
                        //if(Tg[sp].at(current_p)[g]>=bpp.num_base) Tg[sp].at(current_p)[g]=Tg[ss].at(current_s)[g]; ??? 0517

                    }
                }
            }
            
        }        
        ss = sp;
        current_s = current_p;   
    }   
}

/** impute base pair of each internal node on gene tree (except missing nodes) **/
void GTree::Update_Tg(int len, vector<bool> visited, BPP& bpp, bool tosample, vector<int> & nodes, vector<int> & Z, double n_rate, double c_rate, vec log_pi, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv, int start) 
{   
    if(!tosample)
    {
        // 1. sending the lambda msg bottom up & update Tg missing
        for(vector<int>::iterator it = nodes.begin();it<nodes.end();it++)  //
        {
            int ss = *it;
            if(visited[ss] || parent_gene2[ss].size() ==0 ) continue;  // for missing
            //cout<<"inside UpdateTg: node ss="<<ss<<endl;
            double rate = 1;
            if(Z[ss] == 2)
            {
                rate = n_rate;
            }else if(Z[ss] == 1)
            {
                rate = c_rate;
            }            
            for(vector<int>::iterator it2 = temp_coal[ss].begin(); it2 < temp_coal[ss].end(); it2++)
            {
                
                int gs = *it2;
                vector<mat> log_TMs = vector<mat>(2);
                for(int cc=0;cc<2;cc++)
                {
                    int chi = children_gene[gs][cc];
                    double height_child = heights_gene[chi] > bpp.heights[ss] ? heights_gene[chi]:bpp.heights[ss];
                    double lentoCal=heights_gene[gs] - height_child; //length to coal in current species
                    
                    log_TMs[cc] = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv); //ACGT->ACGT transition probability
                }
                for(int g = start; g < len; g++){
                    lambda[ss][gs][g].fill(0); //= zeros<mat>(bpp.num_base, 3);
                    int ct = 0; 
                    for(int cc=0;cc<2;cc++){
                        int chi = children_gene[gs][cc];
                        if(missing_gene[chi] || Tg[ss].at(chi)[g] >= bpp.num_base ){
                            ct ++; 
                            continue;
                        }        
                        //this include: gs ->(single line) kid[gs] -> subtree rooted at kid[sg].            
                        lambda[ss][gs][g].col(cc) = BPP::log_multi(log_TMs[cc], lambda[ss][chi][g].col(2));                        
                    }
                    if(ct==2) Tg[ss][gs][g]=bpp.num_base; 
                    lambda[ss][gs][g].col(2) = lambda[ss][gs][g].col(1) + lambda[ss][gs][g].col(0);
                }                
            }
            
            int sp = bpp.parent[ss];
            if(sp < N)  // send message to parent
            {
                for(map<int, int>::iterator it = parent_gene2[ss].begin(); it!= parent_gene2[ss].end(); it++)
                {
                    if(it->first != it->second) continue; //taken care of in temp_coal;
                    int chi = it->first;
                    double height_child = heights_gene[chi] > bpp.heights[ss] ? heights_gene[chi]:bpp.heights[ss];
                    double lentoCal=bpp.heights[sp] - height_child; 
                    
                    mat log_TM = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv);  
                    
                    // rewrite for speed
                    if(missing_gene[chi])
                    {
                        for(int g = start; g < len; g++)
                        {
                            lambda[sp][chi][g].fill(0); 
                            Tg[sp].at(chi)[g] = Tg[ss].at(chi)[g];
                        }
                    }else{
                        for(int g = start; g < len; g++)
                        {
                            lambda[sp][chi][g].fill(0); 
                            if(Tg[ss][chi][g] >= bpp.num_base)
                            {
                                Tg[sp][chi][g] = Tg[ss][chi][g]; 
                                //above: duplicate. Alr done in InitTg. But InitTg is commented out in Gibbs
                            }else{
                                //P(Y| [sp][chi]position -> [ss][chi]position ->...-> gnode chi creation -> subtree rooted at chi)
                                lambda[sp][chi][g].col(2) = BPP::log_multi(log_TM, lambda[ss][chi][g].col(2)); //chi is a passing gnode, so no left & right kid. Only update col 2.
                            }
                        }
                        
                    }
                }
            }
        }
    }else{
        for(int g = start ;g <len; g++){   
            vec prob = BPP::log_sample(lambda[N-1][root][g].col(2) + log_pi);
            
            //3. sequentially sample the tree top-down
            unsigned int n[bpp.num_base];
            gsl_ran_multinomial(RNG, bpp.num_base, 1,prob.memptr(),n);
            int nn;
            for(nn=0; nn<bpp.num_base;nn++)
            {
                if(n[nn]>0) break;
            }
            
            if(nn>bpp.num_base - 1 || nn<0) {
                cout<<g<<endl;
                throw runtime_error("Sample Tg root error");
                
            }            
            Tg[N-1][root][g] = nn ;
        }
        
        vec log_trans_p = zeros<vec>(bpp.num_base);
        vec trans_p(bpp.num_base);
        for(vector<int>::iterator it = nodes.end()-1; it>= nodes.begin(); it--)
        {
            int ss = *it;
            if(ss < S) continue;
            
            int sp = bpp.parent[ss];
            
            double rate = 1;
            if(Z[ss] == 2)
            {
                rate = n_rate;
            }else if(Z[ss] == 1)
            {
                rate = c_rate;
            }
            
            for(map<int, int>::iterator it = parent_gene2[ss].begin(); it != parent_gene2[ss].end(); it++)
            {
                if(it->first == it->second){ 
                    //second is in species sp. first is in species ss. so here is updating a gnode in ss. gnode can be a new node in ss, or from ss-descendent 
                    int gs = it->first;
                    if(missing_gene[gs]) continue;
                    double height_child = heights_gene[gs] > bpp.heights[ss] ? heights_gene[gs]:bpp.heights[ss]; //Han: whether gs is generated before entering ss, or during ss branch.
                    double lentoCal=bpp.heights[sp] - height_child;
                    if(lentoCal<0){lentoCal=1e-9;}
                    mat log_TM = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv);

                    for(int g = start; g < len; g++)
                    {
                        if(Tg[ss].at(gs)[g] >= bpp.num_base) continue;
                        
                        log_trans_p = log_TM.col(Tg[sp][gs][g]);
                        log_trans_p += lambda[ss][gs][g].col(2); //prior conditional+likelihood
                        
                        trans_p = BPP::log_sample(log_trans_p);
                        
                        unsigned int n[bpp.num_base];
                        gsl_ran_multinomial(RNG, bpp.num_base, 1,trans_p.memptr(),n);
                        
                        int nn;
                        for(nn=0; nn<bpp.num_base;nn++)
                        {
                            if(n[nn]>0) break;
                        }
                        if(nn>bpp.num_base - 1 || nn<0)
                        {
                            throw runtime_error("Sample Tg error");
                        }
                        Tg[ss][gs][g] = nn;
                    }
                    
                }
            }
            //cout<<"end Tg sampling of sp: ss="<<ss<<". temp_coal.size="<<temp_coal[ss].size()<<endl;
           
            if(temp_coal[ss].size() == 0) continue;
            for(vector<int>::iterator it2 = temp_coal[ss].end()-1; it2>= temp_coal[ss].begin(); it2--) 
            {//Han: sample ACGT for internal blue nodes. These are needed, as their ancester will be on sp (red node). So need their info, to sample red node
                int gp = *it2;
                for(int cc =0; cc<2; cc++)
                {
                    int gs = children_gene[gp][cc];
                    if(missing_gene[gs]) continue;
                    double height_child = heights_gene[gs] > bpp.heights[ss] ? heights_gene[gs]:bpp.heights[ss];
                    double lentoCal=heights_gene[gp] - height_child;
                    if(lentoCal<0){lentoCal=1e-9;}
                    mat log_TM = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv);

                    for(int g = start; g < len; g++)
                    {
                        if(Tg[ss].at(gs)[g] >= bpp.num_base) continue;
                        log_trans_p = log_TM.col(Tg[ss][gp][g]);//why can it be 4??? Possibly from Samle_Tree2
                        log_trans_p += lambda[ss][gs][g].col(2);                        
                        trans_p = BPP::log_sample(log_trans_p);
                        
                        unsigned int n[bpp.num_base];
                        gsl_ran_multinomial(RNG, bpp.num_base, 1,trans_p.memptr(),n);
                        
                        
                        int nn;
                        for(nn=0; nn<bpp.num_base;nn++)
                        {
                            if(n[nn]>0) break;
                        }
                        if(nn>bpp.num_base - 1 || nn<0)
                        {
                            throw runtime_error("Sample Tg error");
                        }
                        Tg[ss][gs][g] = nn; //is gs is not a new node. then here is sampling a red node (gnode from below, its seq info at speciation event of ss)
                    }
                }
            }        
        }
    }   
}

void GTree::Simulate_Tg(int len, BPP& bpp, vector<int> & nodes, vector<int> & Z, double n_rate, double c_rate, vec log_pi, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv, int start)
{
    
    // sample root
    for(int g = start ;g <len; g++)
    {
        vec prob = BPP::log_sample(log_pi);
        //cout<<prob.t() <<endl;
        
        //3. sequentially sample the tree top-down
        unsigned int n[bpp.num_base];
        gsl_ran_multinomial(RNG, bpp.num_base, 1,prob.memptr(),n);
        int nn;
        for(nn=0; nn<bpp.num_base;nn++)
        {
            if(n[nn]>0) break;
        }
        
        if(nn>bpp.num_base - 1 || nn<0) {
            cout<<g<<endl;
            throw runtime_error("Sample Tg root error");            
        }        
            Tg[N-1][root][g] = nn ;
    }
    
    
    vec log_trans_p = zeros<vec>(bpp.num_base);
    vec trans_p(bpp.num_base);
    for(vector<int>::iterator it = nodes.end()-1; it>= nodes.begin(); it--)
    {
        int ss = *it;
        
        int sp = bpp.parent[ss];
        
        double rate = 1;
        if(Z[ss] == 2)
        {
            rate = n_rate;
        }else if(Z[ss] == 1)
        {
            rate = c_rate;
        }
        
        
        for(map<int, int>::iterator it = parent_gene2[ss].begin(); it != parent_gene2[ss].end(); it++)
        {
            if(it->first == it->second){
                int gs = it->first;
                if(missing_gene[gs]) continue;
                double height_child = heights_gene[gs] > bpp.heights[ss] ? heights_gene[gs]:bpp.heights[ss];
                double lentoCal=bpp.heights[sp] - height_child;
                if(lentoCal<0){lentoCal=1e-9;}
                mat log_TM = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv);

                for(int g = start; g < len; g++)
                {
                    if(Tg[ss].at(gs)[g] >= bpp.num_base) continue;
                    
                    log_trans_p = log_TM.col(Tg[sp][gs][g]);
                    
                    trans_p = BPP::log_sample(log_trans_p);
                    
                    unsigned int n[bpp.num_base];
                    gsl_ran_multinomial(RNG, bpp.num_base, 1,trans_p.memptr(),n);
                    
                    int nn;
                    for(nn=0; nn<bpp.num_base;nn++)
                    {
                        if(n[nn]>0) break;
                    }
                    if(nn>bpp.num_base - 1 || nn<0)
                    {
                        throw runtime_error("Sample Tg error");
                    }
                    Tg[ss][gs][g] = nn;
                }
                
            }
        }
        
        if(temp_coal[ss].size() ==0) continue;
        
        for(vector<int>::iterator it2 = temp_coal[ss].end()-1; it2>= temp_coal[ss].begin(); it2--)
        {
            int gp = *it2;
            for(int cc =0; cc<2; cc++)
            {
                int gs = children_gene[gp][cc];
                if(missing_gene[gs]) continue;
                double height_child = heights_gene[gs] > bpp.heights[ss] ? heights_gene[gs]:bpp.heights[ss];
                double lentoCal=heights_gene[gp] - height_child;
                if(lentoCal<0){lentoCal=1e-9;}
                mat log_TM = bpp.getlogTMc(lentoCal, rate, c_eigenvec, c_eigenval, c_eigeninv);

                for(int g = start; g < len; g++)
                {
                    if(Tg[ss].at(gs)[g] >= bpp.num_base) continue;
                    
                    log_trans_p = log_TM.col(Tg[ss][gp][g]);
                    
                    trans_p = BPP::log_sample(log_trans_p);
                    
                    unsigned int n[bpp.num_base];
                    gsl_ran_multinomial(RNG, bpp.num_base, 1,trans_p.memptr(),n);
                    
                    
                    int nn;
                    for(nn=0; nn<bpp.num_base;nn++)
                    {
                        if(n[nn]>0) break;
                    }
                    if(nn>bpp.num_base - 1 || nn<0)
                    {
                        throw runtime_error("Sample Tg error");
                    }
                    Tg[ss][gs][g] = nn;
                }
            }
        }        
    }       
}

//Han: sample branch lengths 
bool GTree::Sample_BranchLen(double delta,int gnode, int indicator, BPP & bpp, vector<int> & Z, double n_rate, double c_rate, double consToMis, double nconsToMis, vector<int> lens, vec & loglik,  vec & logp_Z, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv, double & curGTprior, vec & log_pi)
{
    bool braccept=false;
    //double curGTPrior=priorTree(bpp);
    double c0_lens=heights_gene[gnode]-heights_gene[children_gene[gnode][0]]; 
    double c1_lens=heights_gene[gnode]-heights_gene[children_gene[gnode][1]];
    int cid, cnode, sibnode;
    double current_lens;
    //sample the shorter branch, since both kid's branch length depends on the shorter branch. They move exactly the same amount if no topological change.
    if(c0_lens<c1_lens){ 
        //cid=0;
        cnode=children_gene[gnode][0];
        sibnode=children_gene[gnode][1];
        current_lens=c0_lens;
        //sib_lens=c1_lens;
    }else{
        //cid=1;
        cnode=children_gene[gnode][1];
        sibnode=children_gene[gnode][0];
        current_lens=c1_lens;
        //sib_lens=c0_lens;
    }

 
    //Dec2021 Update
    //1. find the max range branch can go:
    //double upperlevel1= parent_gene[gnode]<N ? heights_gene[parent_gene[gnode]]:INFINITY; //max and min value of heights[gnode]
    //double lowerlevel1=heights_gene[cnode]; //just for initialization

    /*
    int ss=gnode;
    while(children_gene[ss][cid]!=-1)
    {
        ss = children_gene[ss][cid];
    }   
    int sp = bpp.parent[ss];
    int ss1=ss; //species for cnode;
    bool tofind1=false;
    if(bpp.heights[sp]<heights_gene[cnode]){
        tofind1=true; //need to find the species that cnode resides in
    }   
    int css=-1; //species in which cnode and sib can coalesce begins. This will take care of the leaf issue.
    map<int, int>::iterator it0=parent_gene2[ss].find(cnode);
    map<int, int>::iterator it1=parent_gene2[ss].find(sibnode);
    if(it0 !=parent_gene2[ss].end() && it1 !=parent_gene2[ss].end()){
        css=ss;
        lowerlevel1=max(bpp.heights[css], heights_gene[cnode]);    
    }
    
    while(sp < N && bpp.heights[sp] <= heights_gene[gnode]) //coalescent strictly higher than speciation
    {   
        
        ss = sp ;        
        sp = bpp.parent[ss];
        if(css==-1){
            it0=parent_gene2[ss].find(cnode);
            it1=parent_gene2[ss].find(sibnode);
            if(it0 !=parent_gene2[ss].end() && it1 !=parent_gene2[ss].end()){
                css=ss;
                lowerlevel1=max(bpp.heights[css], heights_gene[cnode]);
            }
        } 
        if(sp<N){
            if(tofind1 && bpp.heights[sp]>heights_gene[cnode]){
                ss1=ss;
                tofind1=false; 
            }
        }else{
            if(tofind1){
                ss1=ss;
                tofind1=false;
            }
        }             
    } //ss is the species that gnode resides in.
     */
    
    int ss1 = findSp(cnode, bpp); //species for cnode;
    int ss = ss1; //species for gnode;
    int css = ss1; //species for common ancestor of sibnode and cnode;
    int sp = ss; //bpp.parent[ss];
    bool findcss = false;
    // if sp == N, ss1 = ss = N-1;
    while(sp < N  && bpp.heights[sp] < heights_gene[gnode]) // <= ?? 0517
    {
        ss = sp;
        sp = bpp.parent[ss];
        if(!findcss && (parent_gene2[ss].find(sibnode) != parent_gene2[ss].end()))
        {
            css = ss;
            findcss = true;
        }
       
    }
    
    if(!findcss) {
        cout << heights_gene[gnode] <<endl;
        cout << bpp.heights[sp] <<endl;
        cout << bpp.heights[ss] <<endl;
        throw("findcss error");
    }

    
    vector<int>::iterator it = find(temp_coal[ss].begin(), temp_coal[ss].end(), parent_gene[gnode]);
    double upperlevel1;
    //if((it !=temp_coal[ss].end()) || (bpp.move_br[ss]==1)){ 
    if((it !=temp_coal[ss].end())){
        //if parent_gene(gnode) is in the same sp as current pa(cnode), then coal max is pg(gnode) height
        //or if not in, but ss br on sp tree is too short. Max coal is still based on pg(gnode) //17Dec remove this condition. Always within sp move 
        upperlevel1= parent_gene[gnode]<N ? heights_gene[parent_gene[gnode]]:INFINITY; //max and min value of heights[gnode]
    }else{
        upperlevel1 = sp < N ? bpp.heights[sp]:INFINITY;
    }

    if(gnode==root){
        upperlevel1=10*bpp.thetas[N-1]/2; //there is risk upperlevel1 < lowerlevel. 10* is large enough for most cases
        if(upperlevel1<heights_gene[root]) upperlevel1=heights_gene[root]+5*bpp.thetas[N-1]/2; //consistent with Sample_tree2 def of node N height.
    }
    
    double lowerlevel1=max(bpp.heights[css], heights_gene[cnode]);
    
    double upperlevel=min(upperlevel1,heights_gene[gnode]+delta*bpp.thetas[ss]/2);
    double lowerlevel=max(lowerlevel1,heights_gene[gnode]-delta*bpp.thetas[ss]/2);

    double difL=upperlevel-lowerlevel;
    if(difL==0) return(braccept);
    
    
    double prop_lens = max(gsl_rng_uniform(RNG), 1e-5)*difL; //new branch length
    double prop_height = prop_lens + lowerlevel; //proposed height for gnode;
    
    
    //update Tg, lambda: if gnode is above ss, update for cnode and sibnode; otherwise delete
    int ss1pa = css;
    ss1 = css;
    bool findcol = false;
    int newcol = ss1;
    
    mat log_TM;
    do
    {
        ss1 = ss1pa;
        ss1pa=bpp.parent[ss1];
        
        double rate=1;
        if(Z[ss1]==2){
            rate=n_rate;
        }else if(Z[ss1]==1){
            rate=c_rate;
        }
        
        if((ss1pa == N || bpp.heights[ss1pa] >=  prop_height) && !findcol) // get coalescent
        {
            findcol = true;
            newcol = ss1;
            // add gnode
            if(ss != ss1)
            {
                parent_gene2[ss1][cnode] = gnode;
                parent_gene2[ss1][sibnode] = gnode;
                
                if(ss1 < ss)// if new coal higher than old, parent_gene2[ss1][gnode] not changed
                {
                    lambda[ss1].insert(make_pair(gnode,vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
                    //if(ss1pa < N) lambda[ss1pa].insert(make_pair(gnode,vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
                    Tg[ss1].insert(make_pair(gnode,vector<int>(GG,-2)));
                    parent_gene2[ss1][gnode] = gnode;
                }else{
                    
                    //lambda[ss].erase(gnode);
                    //Tg[ss].erase(gnode);
                   // parent_gene2[ss].erase(gnode);
                    
                }
            }
            
            heights_gene[gnode]=prop_height;
            //vector<int>::iterator it = find(temp_coal[ss].begin(), temp_coal[ss].end(), gnode);
            //if(it != temp_coal[ss].end())
            temp_coal[ss].erase(std::remove(temp_coal[ss].begin(), temp_coal[ss].end(), gnode), temp_coal[ss].end());
            
            vector<int>::iterator tempit = temp_coal[ss1].begin();
            for(; tempit != temp_coal[ss1].end(); tempit++)
            {
                if(heights_gene[*tempit] > prop_height)
                {
                    temp_coal[ss1].insert(tempit,gnode);// Han: insert branchp. temp_coal: small height front.
                    break;
                }
            }
            if(tempit == temp_coal[ss1].end()) temp_coal[ss1].push_back(gnode);
            
            for(int i = 0; i<2; i++)
            {
                int cc = children_gene[gnode][i];
                double height_part=heights_gene[cc]>bpp.heights[ss1]?heights_gene[cc]:bpp.heights[ss1];
                double lentoCal = max(prop_height-height_part, 1e-9) ;
                log_TM=bpp.getlogTMc(lentoCal,rate,c_eigenvec,c_eigenval,c_eigeninv);
                
                for(int g=0; g<GG; g++){
                    if(Tg[ss1][cc][g] < bpp.num_base) {
                        lambda[ss1][gnode][g].col(i) = BPP::log_multi(log_TM, lambda[ss1].at(cc)[g].col(2));
                    }
                }
                
            }
            
//            if(ss1pa < N)
//            {
//                double lentoCal = max(bpp.heights[ss1pa] - heights_gene[gnode], 1e-9) ;
//                log_TM=bpp.getlogTMc(lentoCal,rate,c_eigenvec,c_eigenval,c_eigeninv);
//            }
            
            for(int g=0; g<GG; g++){
                int ct = 0;
                for(int i = 0; i<2; i++)
                {
                    int cc = children_gene[gnode][i];
                    if(missing_gene[cc] || Tg[ss1][cc][g] >= bpp.num_base) {
                        ct ++;
                    }
                }
                if(ct==2) Tg[ss1][gnode][g] = bpp.num_base;
                else{
                    lambda[ss1][gnode][g].col(2) = lambda[ss1][gnode][g].col(1) + lambda[ss1][gnode][g].col(0);
                    //if(ss1pa < N) lambda[ss1pa][gnode][g].col(2) = BPP::log_multi(log_TM, lambda[ss1].at(gnode)[g].col(2));
                }

            }
            
        }else if(findcol && ss1 <= ss) // new coal is lower than ss: remove the branches before ss and above the new coalescent node
        {
            //if(ss1 < ss)
            //{
                for(int i = 0; i<2; i++)
                {
                    int cc = children_gene[gnode][i];
                    lambda[ss1].erase(cc);
                    Tg[ss1].erase(cc);
                    parent_gene2[ss1].erase(cc);
                }
            //}

            if(ss1 < ss)
            {
                lambda[ss1].insert(make_pair(gnode,vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
                Tg[ss1].insert(make_pair(gnode,vector<int>(GG,-2)));
                parent_gene2[ss1][gnode] = gnode; //else parent_gene2[ss1pa][gnode] = parent_gene2[ss][gnode];
            }
            
//            double lentoCal = max(bpp.heights[ss1pa]- max(bpp.heights[ss1],  heights_gene[gnode]), 1e-9) ;
//            mat log_TM=bpp.getlogTMc(lentoCal,rate,c_eigenvec,c_eigenval,c_eigeninv);
//
//            if(missing_gene[gnode])
//            {
//                std::fill(Tg[ss1pa][gnode].begin(), Tg[ss1pa][gnode].end(), bpp.num_base);
//            }else{
//
//                for(int g=0; g<GG; g++){
//                    if(Tg[ss1][gnode][g] >= bpp.num_base) {
//                        Tg[ss1pa][gnode][g] = bpp.num_base;
//                    }else{
//                        lambda[ss1pa][gnode][g].col(2) = BPP::log_multi(log_TM, lambda[ss1].at(gnode)[g].col(2));
//                    }
//                }
//            }
            
        }else if(!findcol && ss1 >= ss) //new coal is higher than ss: for the part that branch grows but not coalesced yet. //6-May: add "=" ??? why 0517
        {
            
            double height_part=heights_gene[cnode]>bpp.heights[ss1]?heights_gene[cnode]:bpp.heights[ss1];
            double lentoCal = max(bpp.heights[ss1pa]-height_part, 1e-9) ;
            
            mat log_TM=bpp.getlogTMc(lentoCal,rate,c_eigenvec,c_eigenval,c_eigeninv);
            
            for(int i = 0; i<2; i++)
            {
                int cc = children_gene[gnode][i];
                lambda[ss1pa].insert(make_pair(cc,vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
                Tg[ss1pa].insert(make_pair(cc,vector<int>(GG,-2)));
                
                if(missing_gene[cc])
                {
                    std::fill(Tg[ss1pa][cc].begin(), Tg[ss1pa][cc].end(), bpp.num_base);
                }else{
                
                    for(int g=0; g<GG; g++){
                        if(Tg[ss1][cc][g] >= bpp.num_base) {
                            Tg[ss1pa][cc][g] = bpp.num_base;
                        }else{
                            lambda[ss1pa][cc][g].col(2) = BPP::log_multi(log_TM, lambda[ss1].at(cc)[g].col(2));
                        }
                    }
                }
                parent_gene2[ss1][cc] = cc;
            }
            
            parent_gene2[ss1].erase(gnode);
            lambda[ss1].erase(gnode);
            Tg[ss1].erase(gnode);
            
        }else if(findcol && ss1 > ss) break;
            
    }while(ss1pa<N); //now ss1 is the species that cnode and sibnode coalesce.
   
    
    /*
    int start_ss=Remove_branch(cnode, ss1, bpp);
    Update_Lambda(start_ss, sibnode, bpp, Z, n_rate, c_rate,c_eigenvec, c_eigenval, c_eigeninv);
    
    //update Tg, lambda for the part that branch grows but not coalesced yet.
    int ss1_ori=ss1;
    int ss1pa=bpp.parent[ss1];
    while(ss1pa<N && bpp.heights[ss1pa]<= prop_height){//6-May: add "="
        double rate=1;
        if(Z[ss1]==2){
            rate=n_rate;
        }else if(Z[ss1]==1){
            rate=c_rate;
        }
        double height_part=heights_gene[cnode]>bpp.heights[ss1]?heights_gene[cnode]:bpp.heights[ss1];
        double lentoCal=bpp.heights[ss1pa]-height_part;
        if(lentoCal<0){lentoCal=1e-9;}
        mat log_TM=bpp.getlogTMc(lentoCal,rate,c_eigenvec,c_eigenval,c_eigeninv);
        lambda[ss1pa].insert(make_pair(cnode,vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
        Tg[ss1pa].insert(make_pair(cnode,vector<int>(GG,-2)));
        for(int g=0; g<GG; g++){
            if(missing_gene[cnode] || Tg[ss1][cnode][g] >= bpp.num_base) {
                Tg[ss1pa][cnode][g] = bpp.num_base;
            }else{                
                lambda[ss1pa][cnode][g].col(2) = BPP::log_multi(log_TM, lambda[ss1].at(cnode)[g].col(2));
            }
        }
        parent_gene2[ss1][cnode]=cnode;
        
        ss1=ss1pa;
        ss1pa=bpp.parent[ss1];
    } //now ss1 is the species that cnode and sibnode coalesce.
    double rate=1;
    if(Z[ss1]==2){
        rate=n_rate;
    }else if(Z[ss1]==1){
        rate=c_rate;
    }

    heights_gene[gnode]=prop_height;
    Graft(cnode, ss1, sibnode, prop_height, bpp, rate, c_eigenvec, c_eigenval, c_eigeninv);
     */
    
    
    Update_Lambda(newcol, gnode, bpp, Z, n_rate, c_rate, c_eigenvec, c_eigenval, c_eigeninv);
    
    
    // MH move   
    double r = gsl_rng_uniform(RNG);
    vector<double> loglik_new = vector<double>(2,0);
    for(int g = 0; g< lens[indicator]; g++)
    {
        double temp = BPP::log_exp_sum(lambda[N-1][root][g].col(2) + log_pi);
        loglik_new[0] += temp;
        if(g < lens[1]) loglik_new[1] += temp;
    }
    
    double logp_Z_new = get_logpZ1(bpp, Z, consToMis, nconsToMis);

    double cur_den=-log(upperlevel-lowerlevel);
    double lowerlevel2=max(lowerlevel1,prop_height-delta*bpp.thetas[ss]/2);
    double upperlevel2=min(upperlevel1,prop_height+delta*bpp.thetas[ss]/2);
    double pro_den=-log(upperlevel2-lowerlevel2);
    double delta_loglik = loglik_new[indicator] - loglik[indicator]+pro_den-cur_den;
    
    if(logp_Z[indicator]!=0) // for iter =0 , indicator =1
    {
        delta_loglik += logp_Z_new - logp_Z[indicator];
    }
    double propGTprior=priorTree(bpp);
    delta_loglik+=propGTprior-curGTprior;
    
    if(log(r) < delta_loglik)
    {
        braccept=true;
        loglik = loglik_new;
        if(logp_Z[1]==0) // for iter = 0
        {
            logp_Z[0] = logp_Z_new;
        }else{
            logp_Z[0] = logp_Z[1] = logp_Z_new;
        }
        curGTprior = propGTprior;
    }   
    return(braccept);
}


//compares average difference between two Tg sequence.
double GTree::CompareTg(vector<int> Tg1, vector<int> Tg2, BPP& bpp){ //requires Tg1 not all -2/num_base.
    double diff=0;
    int countMissing=GG;
    for(int i=0; i<GG; i++){
        if(Tg2[i]==bpp.num_base || Tg2[i]==-2){
            continue;
        }else{
            if(Tg1[i] != Tg2[i]) diff=diff+1;
            countMissing--;
        }    
    }
    if(countMissing==GG){ //if Tg2 all missing, return -ve.
        diff=-10;
    }else{ 
         diff=(double) diff/(GG-countMissing);
    } 
    return(diff);
}

void GTree::printSptree(BPP& bpp, int l){
    
    if(l==1){
    cout << "print sp.tree" << endl;
    for (int s = 0; s < N; s++){
        cout <<bpp.nodes_names[s] <<": " << s << " -> " <<  bpp.parent[s] << " (" << bpp.children[s][0] << ", " << bpp.children[s][1]<< "), ht=" << bpp.heights[s] << "; ";
        if(s % 3 ==2) cout << endl;
    }
    }

    cout<<"\n\ngenetree root is"<<root<<endl;

    cout<<"\nvar_br_nodes are:";
    for(int i=0; i<var_br_node.size(); i++) cout<<var_br_node[i]<<" ";
    cout<<endl;
    
    cout << "\nprint gt:" << endl;
    for (int s = 0; s < N; s++){
        cout << s << " -> " <<  parent_gene[s] << " (" << children_gene[s][0] << ", " << children_gene[s][1] << "), ht=" << heights_gene[s]  << "; ";
        if(s % 3 ==2) cout << endl;
    }
    
    cout << "\n\nprint missing gt:" << endl;
    for (int s = 0; s < N; s++){
        if(missing_gene[s]) cout << s << "; " ;
    }
    
    cout << "\n\nprint pg2 info:" << endl;
    for (int s = 0; s < N; s++){
        cout << "sp=" << s << ":";
        for (map<int, int>::iterator it = parent_gene2[s].begin(); it != parent_gene2[s].end(); it++){
            cout << it->first << " -> " << it->second << ", ";
        }
        cout << endl;
    }
    
    cout << "\n\nprint temp_coal info:"<<endl;
    for(int s=S; s<N; s++){
        if(temp_coal[s].size()>0){
            cout<<"temp_coal["<<s<<"]=";
            for (vector<int>::iterator it = temp_coal[s].begin(); it != temp_coal[s].end(); it++){
            int temp_g=*it;
             cout<<temp_g<<", ";
            }
           
        }
        if(s % 3 ==2)  cout<<endl;
    }
    cout<<endl;
   /*
    cout<<"print Tg info:"<<endl;
    for(int s=0; s<N; s++){
        cout <<"sp="<< s << ":";
        for(map<int,vector<int>>::iterator it = Tg[s].begin(); it != Tg[s].end(); it++){
            int temp_e=it->first;
            //cout<<"Tg element: "<<temp_e<<", size="<<gtree.Tg[s][temp_e].size()<<endl;
            cout<<"Tg element:"<<temp_e<<": ";            
            for(int i=0; i<l; i++){
                if(i%5==0) cout<<"["<<i<<"]";
                cout<<Tg[s][temp_e][i];
            }
            cout<<endl;
        }
    }
    */
}

bool GTree::Sample_tree2(int branch, int indicator, BPP & bpp, vector<int> & Z, double n_rate, double c_rate, double consToMis, double nconsToMis, vector<int> lens, vec & loglik,  vec & logp_Z, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv, vec& log_pi)
{
    bool topaccept = false;
    //if(missing_gene[branch]){
    //    topaccept=Sample_tree(indicator,bpp, Z, n_rate,c_rate,consToMis, nconsToMis,lens,loglik,logp_Z,c_eigenvec, c_eigenval, c_eigeninv, log_pi);
    //    return(topaccept);
    //} //if too many missing genes, no info available. So just like sample from prior.

    //bool iscoal = false;
    
    double curGprior=priorTree(bpp);
    // get the lineage on species tree
    int ss =findSp(branch, bpp); //Nov21: ss is the species where gene node branch sits. i.e., possible starting species of coalescence
   
    int branchp = parent_gene[branch];
    int branchsib = children_gene[parent_gene[branch]][1 - childID_gene[branch]];
    int branchpp=parent_gene[branchp];
    double curNHeight=10*bpp.thetas[N-1]/2; //now it's a fixed height, fixed length from sptree.root. Potential risk: node exceed this.
    if(curNHeight<heights_gene[root]){
        curNHeight=heights_gene[root]+5*bpp.thetas[N-1]/2;
    }

    //record sib H seq. sib may be a new brpp. but H get removed in Remove_Branch.
    //int sib_sp=findSp(branchsib, bpp);
   
    //vector<int> H_sib=Tg[sib_sp][branchsib];
    //bool miss_brsib=missing_gene[branchsib];

    //get lowheight under current config, for MH cal.
    double curlowheight;
    int ss_p=ss; //ss_p where they meet
    while(ss_p <N){
        map<int, int>:: iterator it=parent_gene2[ss_p].find(branchsib);
        if(it != parent_gene2[ss_p].end()){
            curlowheight=max(bpp.heights[ss_p],max(heights_gene[branch],heights_gene[branchsib]));
            break;
        }
        ss_p=bpp.parent[ss_p];
    } //ss_p is the first specie that branch and its current sib: branchsib can coalesce.
    
    //30Jun: to be consistent with brlen sampling, upper height in root: 10*mean time; if not root: should be grand-pa height.
    double curuppheight=branchpp<N ? heights_gene[branchpp] : curNHeight;
    double cur_unif_len=-log(curuppheight-curlowheight);
    
    int start_ss = Remove_branch(branch, ss, bpp);

    ////////////////////////////NEW/////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    vector<int> coal_gnode; //potential_gr-parent + parent_sp
    vector<int> coal_gnode_sp;
    vector<double> coal_rl; //coal_group_node coresponding distance to target (r*t) + transition prob.
    vector<double> coal_p;
    double cur_len=0; //part of r*t from ss to gnode_sp (for ss_sp, this is 0)
    double temp_max=-INFINITY;
    ss_p=ss;//coal can start to happen at branch's specie node ss.
    int temp_kid_ssp = ss; //to record specie kid of ss_p along branch's lineage. Use its length to determine if stop candidate searching
    //Dec21 Update
    bool tmp_findgrp = true;
    vector<int> tmp_sibID; //sibID and upLim: for cases where grd is in the pa sp, but actual coal only happen below. The sibling is fixed to be the lineage that coming from the same tracing lineage. And coal upper limit is capped at sp level.
    vector<double> tmp_upLim;
    vector<int> tmp_pg2; //tmp storage lineages that entering and going out from temp_kids_sp. They are potential sib.testSRC_d
    while((ss_p < N) && tmp_findgrp){ //any species from ss and onwards can be potential sp for coalescent
        double temp_rate = 1; //rate is per species.
        if (Z[ss_p] == 2){
            temp_rate = n_rate;
        }else if (Z[ss_p] == 1){
            temp_rate = c_rate;
        }
        if(temp_coal[ss_p].size()==0){
            if(bpp.parent[ss_p]==N){
                coal_gnode.push_back(N);
                coal_gnode_sp.push_back(ss_p);
                tmp_sibID.push_back(-1);
                tmp_upLim.push_back(-1);
                //if no temp_coal, then only 1 lineage left to get coalesced, and become root.
                //artificially assign an length (2*mean of coal) to brpp.
                //it's the length from kid_brpp (if in root sp) or root.sp to brpp.
                double temp_l=temp_rate*(curNHeight-max(bpp.heights[ss_p],heights_gene[branch]))+cur_len;
                coal_rl.push_back(temp_l); //30Jun upperbound
                if(coal_rl.size()==1){
                    coal_p.push_back(log(1));
                }else{
                    mat tempTM=bpp.getlogTM_len(temp_l,c_eigenvec, c_eigenval, c_eigeninv);
                    double temp_p=0;

                    for(int g=0; g<GG; g++){
                        if(Tg[ss][branch][g]>=bpp.num_base) continue; //if missing, no need to add likelihood of this bp.
                        vec tempTMcol=tempTM.col(Tg[ss][branch][g]);
                        temp_p+= BPP::log_exp_sum(tempTMcol+ log_pi); //log(arma::sum(exp(tempTMcol+bpp.log_pi)));
                    }
                     
                    coal_p.push_back(temp_p);
                    if(temp_p>temp_max) temp_max=temp_p;
                }
                break;
            }else{
                cur_len+=(bpp.heights[bpp.parent[ss_p]]-max(bpp.heights[ss_p],heights_gene[branch]))*temp_rate;

                tmp_pg2.clear(); //no grd parent in this node, but have lineage that can coalesce with branch. Existing lineage's grdpa in ancestor sp.
                for(map<int,int>::iterator it_pg2=parent_gene2[ss_p].begin(); it_pg2 !=parent_gene2[ss_p].end(); it_pg2++){
                    int ind_out = it_pg2->first;
                    if(ind_out != branch ) tmp_pg2.push_back(ind_out);
                }

                temp_kid_ssp = ss_p;
                ss_p=bpp.parent[ss_p];
                continue;
            }
        }

        //1st: if no constrained sampling, then unconstrained candidate grdpa.
        //2nd: if in branch's species, then all higher gnode in the sp can be grdpa.
        //3rd: if kid species is a short br, then all nodes in current sp can be grdpa.
        //4th: kid sp no matter long or short, if have no lineage, then all nodes in current sp can be grdpa.
        if((!bpp.cons_sample) || (temp_kid_ssp==ss_p) || bpp.move_br[temp_kid_ssp] == 1 || (tmp_pg2.size()==0)){
            for (vector<int>::iterator it = temp_coal[ss_p].begin(); it != temp_coal[ss_p].end(); it++){
                int temp_g = *it;
                if (ss_p == ss && heights_gene[temp_g] <= heights_gene[branch]) continue;
                coal_gnode.push_back(temp_g);
                tmp_sibID.push_back(-1);
                tmp_upLim.push_back(-1);
                coal_gnode_sp.push_back(ss_p);

                double temp_l = temp_rate * (heights_gene[temp_g] - max(bpp.heights[ss_p], heights_gene[branch])) + cur_len;

                coal_rl.push_back(temp_l);
                mat tempTM = bpp.getlogTM_len(temp_l, c_eigenvec, c_eigenval, c_eigeninv);

                double temp_p = 0;

                for (int g = 0; g < GG; g++) {
                    if (Tg[ss][branch][g] >= bpp.num_base) continue;
                    // option 1: using Tg
                    //                  vec tempTMcol=tempTM.col(Tg[ss][branch][g]);
                    // option 2: using lambda at ss/branch, but Tg at ss_p/temp_g: not ideal, after branch removed, Tg missing pattern and base pair can be changed
                    vec tempTMcol = BPP::log_multi(tempTM, lambda[ss][branch][g].col(2)); // for option 2
                    if (missing_gene[temp_g] || Tg[ss_p][temp_g][g] >= bpp.num_base) {
                        // this is not quite right, because if missing at temp_g, it only means that the subtree underneath is missing
                        // but at temp_g, it still can obtain information from other lineages.
                        temp_p += BPP::log_exp_sum(tempTMcol + log_pi); //log(arma::sum(exp(tempTMcol+bpp.log_pi)));
                    } else {
                        temp_p += tempTMcol[Tg[ss_p][temp_g][g]];
                    }
                }
                coal_p.push_back(temp_p);
                if (temp_p > temp_max)
                    temp_max = temp_p;
            }
            tmp_pg2.clear(); //no grd parent in this node, but have lineage that can coalesce with branch. Existing lineage's grdpa in ancestor sp.
            for(map<int,int>::iterator it_pg2=parent_gene2[ss_p].begin(); it_pg2 !=parent_gene2[ss_p].end(); it_pg2++){
                int ind_out = it_pg2->first;
                if(ind_out != it_pg2->second) continue;
                if(ind_out != branch ) tmp_pg2.push_back(ind_out);
            }
        }else{ //for situation that actual coal will happen in temp_kid_sp as it is long. But current grd is in ancestor node.
            int grd_ct=0;
            for(int tmp_lin = 0; tmp_lin<tmp_pg2.size();tmp_lin++){
                map<int,int>:: iterator it_coal = parent_gene2[ss_p].find(tmp_pg2[tmp_lin]);
                if(it_coal==parent_gene2[ss_p].end()) cout<<"cannot find tmp_pg2[tmp_lin]="<<tmp_pg2[tmp_lin]<<endl;
                int temp_g = it_coal->second;
                if(tmp_pg2[tmp_lin]==temp_g) continue; //this possible sib doesnt coal in ss_p
                grd_ct+=1;
                coal_gnode.push_back(temp_g);
                tmp_sibID.push_back(tmp_pg2[tmp_lin]);
                tmp_upLim.push_back(bpp.heights[ss_p]);
                coal_gnode_sp.push_back(ss_p);
                double temp_l = temp_rate * (heights_gene[temp_g] - max(bpp.heights[ss_p], heights_gene[branch])) + cur_len;

                coal_rl.push_back(temp_l);
                mat tempTM = bpp.getlogTM_len(temp_l, c_eigenvec, c_eigenval, c_eigeninv);

                double temp_p = 0;

                for (int g = 0; g < GG; g++) {
                    if (Tg[ss][branch][g] >= bpp.num_base) continue;
                    // option 1: using Tg
                    //                  vec tempTMcol=tempTM.col(Tg[ss][branch][g]);
                    // option 2: using lambda at ss/branch, but Tg at ss_p/temp_g: not ideal, after branch removed, Tg missing pattern and base pair can be changed
                    vec tempTMcol = BPP::log_multi(tempTM, lambda[ss][branch][g].col(2)); // for option 2
                    if (missing_gene[temp_g] || Tg[ss_p][temp_g][g] >= bpp.num_base) {
                        // this is not quite right, because if missing at temp_g, it only means that the subtree underneath is missing
                        // but at temp_g, it still can obtain information from other lineages.
                        temp_p += BPP::log_exp_sum(tempTMcol + log_pi); //log(arma::sum(exp(tempTMcol+bpp.log_pi)));
                    } else {
                        temp_p += tempTMcol[Tg[ss_p][temp_g][g]];
                    }
                }
                coal_p.push_back(temp_p);
                if (temp_p > temp_max)
                    temp_max = temp_p;
            }
            if(grd_ct >0){
                tmp_findgrp=false; //stop at this ss_p
            }else{
                tmp_pg2.clear(); //no grd parent in this node, but have lineage that can coalesce with branch. Existing lineage's grdpa in ancestor sp.
                for(map<int,int>::iterator it_pg2=parent_gene2[ss_p].begin(); it_pg2 !=parent_gene2[ss_p].end(); it_pg2++){
                    int ind_out = it_pg2->first;
                    if(ind_out != it_pg2->second) continue;
                    if(ind_out != branch ) tmp_pg2.push_back(ind_out);
                }                
            }
        }
        
        if((bpp.parent[ss_p]==N)){
            double temp_ht=0;
            if(coal_gnode.size()>0){
                temp_ht=heights_gene[coal_gnode[coal_gnode.size()-1]];
            }else temp_ht=max(bpp.heights[ss_p],heights_gene[branch]);
            
            coal_gnode.push_back(N);
            tmp_sibID.push_back(-1);
            tmp_upLim.push_back(-1);
            coal_gnode_sp.push_back(ss_p);
            double temp_l=cur_len+temp_rate*(curNHeight-temp_ht);
            
            coal_rl.push_back(temp_l);  //30Jun: from branch  - N. 
            mat tempTM=bpp.getlogTM_len(temp_l,c_eigenvec, c_eigenval, c_eigeninv);
            double temp_p=0;
            
            for(int g=0; g<GG; g++){
                if(Tg[ss][branch][g]>=bpp.num_base) continue;
                
                // option 1
                //vec tempTMcol=tempTM.col(Tg[ss][branch][g]);
                //temp_p+=log(sum(exp(tempTMcol+bpp.log_pi)));
                
                //option 2
                vec tempTMcol = BPP::log_multi(tempTM, lambda[ss][branch][g].col(2));
                temp_p+= BPP::log_exp_sum(tempTMcol+log_pi);
                
            }

            coal_p.push_back(temp_p);
            if(temp_p>temp_max) temp_max=temp_p;

        }
        
        if(bpp.parent[ss_p] < N){
            cur_len+=(bpp.heights[bpp.parent[ss_p]]-max(bpp.heights[ss_p],heights_gene[branch]))*temp_rate;
        }

        if(ss_p == bpp.root_ingrp){
            if(coal_gnode.size()>0){
                tmp_findgrp=false;
            }
        } 

        temp_kid_ssp = ss_p;
        ss_p=bpp.parent[ss_p];
    } //get all potential brpp for target.

    int get_p;//coal_group_node[0][get_p] is the new branchpp.
    vec sample_p = zeros<vec>(coal_p.size());
    if(coal_gnode.size()==1){
        get_p=0;
        sample_p[0]=1;
    }else{
        sample_p = BPP::log_sample_norm(coal_p);

        unsigned int new_p[coal_gnode.size()];
        gsl_ran_multinomial(RNG, coal_gnode.size(), 1, sample_p.memptr(), new_p);
        for (get_p = 0; get_p < coal_gnode.size(); get_p++){
            if (new_p[get_p] > 0)
                break;
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////END NEW///////////////////////////////////////////////////////
    Update_Lambda(start_ss, branchsib, bpp, Z, n_rate, c_rate,c_eigenvec, c_eigenval, c_eigeninv);     
    ////////////////////////////////NEW//////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    int new_brpp=coal_gnode[get_p];//gnode
    int new_sib;
    double uppheight;
    if(new_brpp==N){
        new_sib=root;
        uppheight=10*bpp.thetas[N-1]/2 > heights_gene[new_sib] ? 10*bpp.thetas[N-1]/2 : heights_gene[new_sib]+5*bpp.thetas[N-1]/2;
    }else if(tmp_sibID[get_p]!=(-1)){ //if grp node is in the pa species of the long specie br where coal should have happended, then must get sibling node from the below long br, but not from the other lineage
        new_sib = tmp_sibID[get_p]; //this gnode is along gnode's lineage, while the othr kid of coal_gnode will be from another lineage
        uppheight=tmp_upLim[get_p];
    }else{
        unsigned long a = gsl_rng_uniform_int(RNG, 2);
        new_sib=children_gene[new_brpp][a]; //gnode
        uppheight=heights_gene[new_brpp];
    }
    //find the species that branch and new_sib can start coalesce. This is the lower bound for brp position
    double lowheight;
    ss_p=ss;
    while(ss_p != bpp.parent[coal_gnode_sp[get_p]]){
        map<int, int>:: iterator it=parent_gene2[ss_p].find(new_sib);
        if(it !=parent_gene2[ss_p].end()){
            lowheight=max(heights_gene[branch],max(heights_gene[new_sib],bpp.heights[ss_p]));
            break;
        }else{
            ss_p=bpp.parent[ss_p];
        }
    }
    double new_unif_len;
    double coal = gsl_rng_uniform(RNG);
    if (coal < 1e-5) coal = 1e-5;
    coal = coal * (uppheight - lowheight) + lowheight;
    heights_gene[branchp] = coal;
    new_unif_len = -log(uppheight - lowheight);
    ////////////////////////////////END Sampling sib//////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
 
    //update Tg, lambda for the part that branch grows but not coalesced yet.
    ss_p=bpp.parent[ss];
    while(ss_p<N && (bpp.heights[ss_p]< heights_gene[branchp])){ // < ?? 0517
        double temp_rate=1;
        if(Z[ss]==2){
            temp_rate=n_rate;
        }else if(Z[ss]==1){
            temp_rate=c_rate;
        }
        double height_part=heights_gene[branch]>bpp.heights[ss]?heights_gene[branch]:bpp.heights[ss];
        double lentoCal=bpp.heights[ss_p]-height_part;
        mat log_TM=bpp.getlogTMc(lentoCal,temp_rate,c_eigenvec,c_eigenval,c_eigeninv);
        lambda[ss_p].insert(make_pair(branch,vector<mat>(GG, zeros<mat>(bpp.num_base, 3))));
        Tg[ss_p].insert(make_pair(branch,vector<int>(GG,-2)));
        for(int g=0; g<GG; g++){
            if(missing_gene[branch] || Tg[ss][branch][g] >= bpp.num_base) {
                Tg[ss_p][branch][g] = bpp.num_base; //should I sample from pi? 
            }else{                
                lambda[ss_p][branch][g].col(2) = BPP::log_multi(log_TM, lambda[ss].at(branch)[g].col(2));
            }
        }
        parent_gene2[ss][branch]=branch;       
        ss=ss_p;
        ss_p=bpp.parent[ss];
    } //now ss is the species that cnode and sibnode coalesce.

    double temp_rate=1;
    if(Z[ss]==2){
        temp_rate=n_rate;
    }else if(Z[ss]==1){
        temp_rate=c_rate;
    }
    Graft(branch, ss,new_sib,heights_gene[branchp], bpp, temp_rate, c_eigenvec, c_eigenval, c_eigeninv);
    Update_Lambda(ss, branchp, bpp, Z, n_rate, c_rate, c_eigenvec, c_eigenval, c_eigeninv);
    // MH move
    double newGprior=priorTree(bpp);
    double r = gsl_rng_uniform(RNG);
    vector<double> loglik_new = vector<double>(2,0);
    for(int g = 0; g< lens[indicator]; g++)
    {
        double temp = BPP::log_exp_sum(lambda[N-1][root][g].col(2) + log_pi);
        loglik_new[0] += temp;
        if(g < lens[1]) loglik_new[1] += temp;
    }
    
    double logp_Z_new = get_logpZ1(bpp, Z, consToMis, nconsToMis);
    double delta_loglik = loglik_new[indicator] - loglik[indicator];
    if(logp_Z[indicator]!=0) // for iter =0 , indicator =1
    {
        delta_loglik += logp_Z_new - logp_Z[indicator];
    }
    delta_loglik+=newGprior-curGprior;
    delta_loglik+=cur_unif_len-new_unif_len;
    //get multinom prob for cur brpp:
    vector<int>:: iterator it_curpp;
    it_curpp=find(coal_gnode.begin(),coal_gnode.end(),branchpp);
    int it_pos=it_curpp-coal_gnode.begin();
    //double curmultiP=coal_gnode_sp[it_pos]; //30Jun
    delta_loglik+=log(sample_p[it_pos])-log(sample_p[get_p]);
    
    if(log(r) < delta_loglik) // accept
    {
        topaccept = true;
        loglik = loglik_new;
        if(logp_Z[1]==0) // for iter = 0
        {
            logp_Z[0] = logp_Z_new;
        }else{
            logp_Z[0] = logp_Z[1] = logp_Z_new;
        }
    }
    return(topaccept);
}


int GTree::findSp(int gnode, BPP& bpp){
    int gnode_sp = gnode;
    while(children_gene[gnode_sp][0]!=-1 || children_gene[gnode_sp][1] != -1){//8Nov, as can be -1 in one kid from out group
        gnode_sp = children_gene[gnode_sp][0];
        if(gnode_sp==-1) gnode_sp = children_gene[gnode_sp][1]; 
    }
    //to add in code to handle internal nodes along outgrp lineage that is missing both kids 
    //if(gnode_sp>= S) //may use parent_gene2 / temp_coal to find
    int ss_pa = bpp.parent[gnode_sp]; //this works if gnode is not among out-group lineages.  //can make restriction to never sample along out-group lineages
    double current_height = heights_gene[gnode];   
    while(ss_pa < N && bpp.heights[ss_pa] < current_height){ // ? '<' 0517
        gnode_sp = ss_pa ;   
        ss_pa = bpp.parent[gnode_sp];
    } 
    return(gnode_sp);
}
