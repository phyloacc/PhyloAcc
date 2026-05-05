#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <algorithm>
#include <vector>
#include <string>
#include <set>
#include <numeric>

using namespace std;

namespace ctnutils
{
    template <typename T>
    vector<size_t> sort_indexes(const vector<T> &v) {
        vector<size_t> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        sort(idx.begin(), idx.end(),
             [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
        return idx;
    }

template <typename T>
struct IdxCompare
{
    IdxCompare(const T target) : target(target) {}
    bool operator()(const size_t a, const size_t b) const
    {
        return target[a].size() > target[b].size();
    }
    const T target;
};

template <typename T>
vector<size_t> sort_indexes(const T & target)
{
    vector<size_t> idx(target.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

    sort(idx.begin(), idx.end(), IdxCompare<T>(target));

    return idx;
}

template<class T>
vector<T> merge(const vector<T> & v1, const vector<T> & v2)
{
    vector<T> v3;
    v3.reserve(v1.size() + v2.size());
    v3.insert(v3.end(), v1.begin(), v1.end());
    v3.insert(v3.end(), v2.begin(), v2.end());
    return v3;
}

template<class T>
vector<T> setdiff(const vector<T> & a, const vector<T> & b)
{
    set<T> diff(a.begin(), a.end());
    for (const T & item : b)
        diff.erase(item);
    return vector<T>(diff.begin(), diff.end());
}

template<class T>
vector<T> intersect(const vector<T> & a, const vector<T> & b)
{
    vector<T> out;
    set<T> lookup(b.begin(), b.end());
    for (const T & item : a)
        if (lookup.count(item))
            out.push_back(item);
    return out;
}

template<class T>
bool isin(const T& x, const vector<T> & v)
{
    return find(v.begin(), v.end(), x) != v.end();
}
}

namespace strutils
{
vector<string> split(string s,char splitchar=' ');
string trim(string str, string chars = " ");
string ToLowerCase(string str);
bool is_number(const std::string& s);
}

#endif
