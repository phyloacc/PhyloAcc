#include "profile.h"
#include "utils.h"

#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;

PhyloProf LoadPhyloProfiles(string profile_path, string segment_path, string segment_ID)
{
    PhyloProf prof;
    prof.G = 0;
    string linestr;

    ifstream in_prof(profile_path.c_str());

    if (!in_prof)
    {
        cerr << "(Error. Cannot open the phylogenetic profile input file: " << profile_path << ")" << endl;
        exit(1);
    }

    string wholeline = "";
    while (!in_prof.eof())
    {
        std::getline(in_prof, linestr);
        linestr = strutils::trim(linestr);
        if (!strncmp(linestr.c_str(), ">", 1))
        {
            string tmp = strutils::trim(linestr.substr(1));
            prof.species_names.push_back(tmp);

            if (prof.G == 0)
                prof.G = wholeline.length();
            else
                assert(wholeline.length() == prof.G);

            if (prof.G > 0)
            {
                wholeline = strutils::ToLowerCase(wholeline);
                prof.X.push_back(wholeline);
            }
            wholeline = "";
        }
        else
        {
            wholeline += strutils::trim(linestr);
        }
    }
    if (prof.G == 0)
        prof.G = wholeline.length();
    wholeline = strutils::ToLowerCase(wholeline);
    prof.X.push_back(wholeline);

    prof.S = prof.species_names.size();

    ifstream in_segment(segment_path.c_str());

    if (!in_segment)
    {
        cerr << "(Error. Cannot open the segment input file: " << segment_path << ")" << endl;
        exit(1);
    }

    while (!in_segment.eof())
    {
        std::getline(in_segment, linestr);
        linestr = strutils::trim(linestr);
        if (linestr == "")
            continue;
        vector<string> line_splits = strutils::split(linestr, '\t');
        if (line_splits.size() < 3)
            break;
        prof.element_names.push_back(line_splits[0]);
        double* tmp = new double[3];
        tmp[0] = atoi(line_splits[1].c_str());
        tmp[1] = atoi(line_splits[2].c_str());
        prof.element_pos.push_back(tmp);
        if (line_splits.size() >= 8)
            prof.element_tree.push_back(line_splits[7]);
    }
    prof.C = prof.element_names.size();

    in_segment.close();

    if (segment_ID != "")
    {
        string segment_path2 = segment_ID + ".txt";
        in_segment.open(segment_path2.c_str());

        if (!in_segment)
        {
            cerr << "(Error. Cannot open the segment input txt file: " << segment_path2 << ")" << endl;
            exit(1);
        }

        while (!in_segment.eof())
        {
            std::getline(in_segment, linestr);
            linestr = strutils::trim(linestr);
            if (linestr == "")
                continue;
            prof.element_id.push_back(linestr);
        }
    }

    in_segment.close();

    return prof;
}
