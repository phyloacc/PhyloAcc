#include "utils.h"

#include <sys/types.h>
#include <dirent.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;

namespace ctnutils
{

vector<string> GenGrayArray(int n)
{
    // base case
    if (n <= 0)
        return vector<string>();

    // 'arr' will store all generated codes
    vector<string> arr;

    // start with one-bit pattern
    arr.push_back("0");
    arr.push_back("1");

    // Every iteration of this loop generates 2*i codes from previously
    // generated i codes.
    int i, j;
    for (i = 2; i < (1<<n); i = i<<1)
    {
        // Enter the prviously generated codes again in arr[] in reverse
        // order. Nor arr[] has double number of codes.
        for (j = i-1 ; j >= 0 ; j--)
            arr.push_back(arr[j]);

        // append 0 to the first half
        for (j = 0 ; j < i ; j++)
            arr[j] = "0" + arr[j];

        // append 1 to the second half
        for (j = i ; j < 2*i ; j++)
            arr[j] = "1" + arr[j];
    }

//    // print contents of arr[]
//    for (i = 0 ; i < arr.size() ; i++ )
//        cout << arr[i] << endl;

    return arr;
}

}

namespace strutils
{

string LoadFileToString(ifstream & in)
{
    istreambuf_iterator<char> beg(in), end;
    string str(beg, end);
    return str;
}

string LoadFileToString(const char *file_path)
{
    ifstream in(file_path);
    istreambuf_iterator<char> beg(in), end;
    string str(beg, end);
    return str;
}

string LoadFileToString(string file_path)
{
    ifstream in(file_path.c_str());
    istreambuf_iterator<char> beg(in), end;
    string str(beg, end);
    return str;
}

vector<string> &split(const string &s, char delim, vector<string> &elems)
{
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim))
    {
//        if (item.length() > 0)      // skip empty token
            elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim)
{
    vector<string> elems;
    return split(s, delim, elems);
}

string ToUpperCase(string str)
{
    transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

string ToLowerCase(string str)
{
    transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

string replace(string value, string const & search, string const & replace)
{
    string::size_type next;

    for(next = value.find(search);next != std::string::npos;next = value.find(search,next))
    {
        value.replace(next,search.length(),replace);   // Do the replacement.
        next += replace.length();                                      // the next search from.
    }

    return value;
}

bool endsWith (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

string itoa(int value)
{
    std::string buf;

    enum { kMaxDigits = 35 };
    buf.reserve( kMaxDigits ); // Pre-allocate enough space.

    int quotient = value;

    // Translating number to string with base:
    do {
        buf += "0123456789abcdef"[ std::abs( quotient % 10 ) ];
        quotient /= 10;
    } while ( quotient );

    // Append the negative sign
    if ( value < 0) buf += '-';

    std::reverse( buf.begin(), buf.end() );
    return buf;
}

}

namespace numutils
{
    int CountDigits(int num)
    {
        if (num==0)
            return 1;
        int num_digits = (num<0)?1:0;
        while(num!=0)
        {
            num_digits ++;
            num /= 10;
        }
        return num_digits;
    }

    string EnoughSpaces(int MAX, int num)
    {
        string spaces = "";
        int num_spaces = MAX-CountDigits(num);
        for(int i=0; i<num_spaces; i++)
            spaces = spaces + " ";
        return spaces;
    }

    bool is_number(const std::string& s)
    {
        std::string::const_iterator it = s.begin();
        while (it != s.end() && std::isdigit(*it)) ++it;
        return !s.empty() && it == s.end();
    }
}

namespace oututils
{
    void RepeatOutputs(int num_spaces, ostream & out, string sp)
    {
        for(int i=0; i<num_spaces; i++)
            out << sp;
    }
}

namespace fileutils
{
    vector<string> getdir(string dir)
    {
        vector<string> files;
        DIR *dp;
        struct dirent *dirp;
        if((dp  = opendir(dir.c_str())) == NULL) {
            cout << "Error opening " << dir << endl;
            exit(1);
        }

        while ((dirp = readdir(dp)) != NULL) {
            files.push_back(dir+"/"+string(dirp->d_name));
        }
        closedir(dp);
        return files;
    }

    bool file_exists(const std::string& name)
    {
        ifstream f(name.c_str());
        if (f.good()) {
            f.close();
            return true;
        } else {
            f.close();
            return false;
        }
    }
}
