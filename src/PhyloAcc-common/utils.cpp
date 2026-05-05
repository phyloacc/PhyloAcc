#include "utils.h"

#include <cctype>
#include <sstream>

using namespace std;

namespace strutils
{
vector<string> split(string s, char splitchar)
{
    vector<string> items;
    string item = "";
    for (unsigned i = 0; i < s.length(); i++)
    {
        if (s[i] == splitchar)
        {
            items.push_back(item);
            item = "";
        }
        else
        {
            item += s[i];
        }
    }
    items.push_back(item);
    return items;
}

string trim(string str, string chars)
{
    str.erase(0, str.find_first_not_of(chars));
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}

string ToLowerCase(string str)
{
    for (unsigned i = 0; i < str.length(); i++)
    {
        str[i] = tolower(str[i]);
    }
    return str;
}

bool is_number(const std::string& s)
{
    std::istringstream iss(s);
    double d;
    return (iss >> d) && (iss >> std::ws).eof();
}
}
