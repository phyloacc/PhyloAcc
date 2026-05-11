#include "profile.h"
#include "utils.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;

namespace {

void ValidateProfileSequenceLength(const string& profile_path,
                                   const string& species_name,
                                   size_t observed_length,
                                   unsigned expected_length)
{
    if (observed_length == expected_length)
        return;

    cerr << "(Error. Sequence length mismatch in phylogenetic profile input file: "
         << profile_path << ". Species " << species_name << " has "
         << observed_length << " sites; expected " << expected_length << ".)" << endl;
    exit(1);
}

}  // namespace

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
            if (!prof.species_names.empty())
            {
                if (prof.G == 0)
                    prof.G = wholeline.length();
                else
                    ValidateProfileSequenceLength(profile_path, prof.species_names.back(),
                                                  wholeline.length(), prof.G);

                if (prof.G > 0)
                {
                    wholeline = strutils::ToLowerCase(wholeline);
                    prof.X.push_back(wholeline);
                }
            }
            string tmp = strutils::trim(linestr.substr(1));
            prof.species_names.push_back(tmp);
            wholeline = "";
        }
        else
        {
            wholeline += strutils::trim(linestr);
        }
    }
    if (!prof.species_names.empty())
    {
        if (prof.G == 0)
            prof.G = wholeline.length();
        else
            ValidateProfileSequenceLength(profile_path, prof.species_names.back(),
                                          wholeline.length(), prof.G);
        wholeline = strutils::ToLowerCase(wholeline);
        prof.X.push_back(wholeline);
    }
    else
    {
        if (prof.G == 0)
            prof.G = wholeline.length();
        wholeline = strutils::ToLowerCase(wholeline);
        prof.X.push_back(wholeline);
    }

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
