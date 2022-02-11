#include "profile.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <cassert>

#include "utils.h"

using namespace std;

// load the phylogenetic profile
PhyloProf LoadPhyloProfiles(string profile_path, string segment_path, string segment_ID)
{
    PhyloProf prof;
    vector<string> seqs;
    prof.G=0;
    string linestr;
    ifstream in_prof(profile_path.c_str());

    if (!in_prof)
    {
        cerr << "(Error. Cannot open the phylogenetic profile input file: " << profile_path << ")" << endl;
        exit(1);
    }


    

    // count the num of species, base pairs and load the profiles
    string wholeline="";
    while(!in_prof.eof())
    {
        std::getline(in_prof, linestr);
        linestr = strutils::trim(linestr);
        if(!strncmp(linestr.c_str(),">", 1)) {
            string tmp = strutils::trim(linestr.substr(1));
            prof.species_names.push_back(tmp);
            
            if(prof.G==0)   prof.G = wholeline.length();
            else assert(wholeline.length() == prof.G);
            
            if(prof.G>0) {
                wholeline =strutils::ToLowerCase(wholeline);
                prof.X.push_back(wholeline);
            }
            wholeline = "";
        }
        else {
            
            wholeline += linestr;
        }
    }
    
    if(prof.G==0)   prof.G = wholeline.length();
    wholeline =strutils::ToLowerCase(wholeline);
    prof.X.push_back(wholeline);
    
    //cout <<prof.G;
    
    prof.S = prof.species_names.size();
  
    
    
     //read in segment size and specific scaling factor
    
    ifstream in_segment(segment_path.c_str());
    
    if (!in_segment)
    {
        cerr << "(Error. Cannot open the segment input file: " << segment_path << ")" << endl;
        exit(1);
    }

    while(!in_segment.eof())
    {
        std::getline(in_segment, linestr);
        linestr = strutils::trim(linestr);
        vector<string> line_splits = strutils::split(linestr, '\t');
        if(line_splits.size()<3) break;
        prof.element_names.push_back(line_splits[0]);
        double* tmp = new double[3];
        tmp[0] = atoi(line_splits[1].c_str());
        tmp[1] = atoi(line_splits[2].c_str());
        //tmp[2] = atof(line_splits[4].c_str());  //add null scale!!
        prof.element_pos.push_back(tmp);
    }
    prof.C = prof.element_names.size();
    
    in_segment.close();
    
    
    // read in ID
    if(segment_ID!="")
    {
    string segment_path2 = segment_ID + ".txt";
    in_segment.open(segment_path2.c_str());
    
    if (!in_segment)
    {
        cerr << "(Error. Cannot open the segment input txt file: " << segment_path2 << ")" << endl;
        exit(1);
    }
    
    while(!in_segment.eof())
    {
        std::getline(in_segment, linestr);
        linestr = strutils::trim(linestr);
        if(linestr=="") continue;
        //vector<string> line_splits = strutils::split(linestr, '\t');
        //if(line_splits.size()<3) break;
        prof.element_id.push_back(linestr);
    }
    }
    
    
    in_segment.close();
   

    return prof;

}


