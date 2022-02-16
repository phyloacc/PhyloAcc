#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED


#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <set>

using namespace std;

namespace ctnutils
{

//template<class T>
//vector< vector<T> > sortbySizeDec(vector< vector<T> > v)
//{
//    for(int i=0; i<v.size(); i++)
//    {
//        for(int j=i+1; j<v.size(); j++)
//        {
//            if (v[i].size() < v[j].size())
//            {
//
//            }
//        }
//    }
//}
    
template <typename T>
struct IdxCompare
{
    const std::vector<T>& target;
        
    IdxCompare(const std::vector<T>& target): target(target) {}
        
    bool operator()(size_t a, size_t b) const { return target[a].size() > target[b].size(); }
};

template <class T>
vector<size_t> sort_indexes(const vector<T> &v, int K) {
   
    
    if(K<0) K =v.size();
    // initialize original index locations
    vector<size_t> idx(K);
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), IdxCompare<T>(v));
    
    return idx;
}
    

//    template <typename T>
//    vector<size_t> sort_indexes(const vector<T> &v, int K) {
//    
//        if(K<0) K =v.size();
//        // initialize original index locations
//        vector<size_t> idx(K);
//        for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
//    
//        // sort indexes based on comparing values in v
//        sort(idx.begin(), idx.end(),
//             [&v](size_t i1, size_t i2) {return v[i1].size() > v[i2].size();});
//        
//        return idx;
//    }



template<class T>
vector<T> merge(const vector<T> & v1, const vector<T> & v2)
{
    vector<T> c;
    c.insert(c.end(), v1.begin(), v1.end());
    c.insert(c.end(), v2.begin(), v2.end());
    return c;
}

template<class T>
vector<T> merge(const vector<T> & v1, const T & v2)
{
    vector<T> c;
    c.insert(c.end(), v1.begin(), v1.end());
    c.push_back(v2);
    return c;
}

vector<string> GenGrayArray(int n);

template<class T>
vector< vector<T> > subsets(const vector<T> & s)
{
    int len = s.size();
    vector<string> arr = ctnutils::GenGrayArray(len);

    vector< vector<T> > results;
    for(int i=0; i<arr.size(); i++)
    {
        vector<T> result;
        for(int j=0; j<len; j++)
        {
            if (arr[i][j] == '1')
            {
                result.push_back(s[j]);
            }
        }
        results.push_back(result);
    }

    return results;
}

// C = A - B
template<class T>
vector<T> SetDiff(const vector<T> & A, const vector<T> & B)
{
    vector<T> C = A;
    for(int i=0; i<C.size(); i++)
    {
        bool find = false;
        for(int j=0; j<B.size(); j++)
            if (C[i] == B[j])
            {
                find = true;
                break;
            }
        if (find)
        {
            C.erase(C.begin()+i);
            i--;
        }
    }
    return C;
}

template<class T>
void DispVector(const vector<T> & vec, string sep = "\t")
{
    if (vec.size()==0)
    {
        cout << "(empty vector)";
        return;
    }

    cout << "(";
    for(unsigned i=0; i<vec.size()-1; i++)
        cout << vec[i] << sep;
    cout << vec[vec.size()-1] << ")";
}

template<class T>
void DispVector(const vector<int> & vec, const vector<T> & ref, string sep = "\t")
{
    if (vec.size()==0)
    {
        cout << "(empty vector)";
        return;
    }

    cout << "(";
    for(unsigned i=0; i<vec.size()-1; i++)
        cout << ref[vec[i]] << sep;
    cout << ref[vec[vec.size()-1]] << ")";
}

template<class T>
void DispVector(T *vec, int size, string sep = "\t")
{
    if (size==0)
    {
        cout << "(empty vector)";
        return;
    }

    cout << "(";
    for(int i=0; i<size; i++)
        cout << vec[i] << sep;
    cout << vec[size-1] << ")";
}

template<class T>
void DispMatrix(vector<vector<T> > & matrix)
{
    for(unsigned i=0; i<matrix.size(); i++)
    {
        if (matrix[i].size()==0)
            cout << "(empty row)";

        for(unsigned j=0; j<matrix[i].size(); j++)
            cout << matrix[i][j] << "\t";

        cout << endl;
    }
}

template<class T>
void DispMatrix(T **matrix, int num_row, int num_col)
{
    for(int i=0; i<num_row; i++)
    {
        for(int j=0; j<num_col; j++)
            cout << matrix[i][j] << "\t";

        cout << endl;
    }
}

template<class T>
void FillVector(T *vec, int n, T elem)
{
	for(int i=0;i<n;i++)
		vec[i] = elem;
}

template<class T>
void FillVector(vector<T> & vec, int n, T elem)
{
	for(int i=0;i<n;i++)
		vec[i] = elem;
}

template<class T>
bool IsVectorSame(const vector<T> & v1, const vector<T> & v2)
{
    if (v1.size()!=v2.size())
        return false;
    for(unsigned i=0; i<v1.size(); i++)
        if (v1[i]!=v2[i])
            return false;
    return true;
}

template<class T>
bool IsVectorSame(T *v1, T *v2, int n)
{
    for(int i=0; i<n; i++)
        if (v1[i]!=v2[i])
            return false;
    return true;
}

template<class T>
bool Contains(const vector<T> & vec, T elem)
{
    for(int i=0; i<vec.size(); i++)
    {
        if (vec[i]==elem)
            return true;
    }
    return false;
}


template<class T>
bool Contains(const vector<T> & vec, const vector<T> & vec2)
{
    for(int i=0; i<vec2.size(); i++)
        if (!Contains(vec, vec2[i]))
            return false;
    return true;
}


template<class T>
struct FindReturn
{
    int  pos;
    T    val;
};

template<class T>
FindReturn<T> Find(const vector<T> & vec, T t)
{
    FindReturn<T> fr;
    fr.pos = -1;
    for(int i=0; i<vec.size(); i++)
    {
        if (vec[i] == t)
        {
            fr.pos = i;
            break;
        }
    }
    return fr;
}

template<class T>
FindReturn<T> FindMax(const vector<T> & vec, T min)
{
    FindReturn<T> fr;
    fr.pos = -1;
    for(int i=0; i<vec.size(); i++)
    {
        if (vec[i] > min)
        {
            fr.pos = i;
            fr.val = vec[i];
            min    = vec[i];
        }
    }
    return fr;
}

template<class T>
FindReturn<T> FindMin(const vector<T> & vec, T max)
{
    FindReturn<T> fr;
    fr.pos = -1;
    for(int i=0; i<vec.size(); i++)
    {
        if (vec[i] < max)
        {
            fr.pos = i;
            fr.val = vec[i];
            max    = vec[i];
        }
    }
    return fr;
}

}

namespace strutils
{

string LoadFileToString(const char *file_path);
string LoadFileToString(string file_path);
string LoadFileToString(ifstream & in);

inline string trim(string str, string ts = " \t\n")
{
    str.erase(0, str.find_first_not_of(ts));       //prefixing spaces
    str.erase(str.find_last_not_of(ts)+1);         //surfixing spaces

    return str;
}

vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string>  split(const string &s, char delim);
    
    
string ToUpperCase(string str);

string ToLowerCase(string str);

string replace(string value, string const & search, string const & replace);

bool endsWith (std::string const &fullString, std::string const &ending);

string itoa(int value);
}

namespace numutils
{

int CountDigits(int num);

string EnoughSpaces(int MAX, int num);

bool is_number(const std::string& s);

}

namespace oututils
{

void RepeatOutputs(int num_spaces, ostream & out = cout, string sp = " ");

}


namespace statutils
{
// adjusted rand index
double ari(vector<int> I1, vector<int> I2);

double ari(vector<vector<int> > I1, vector<vector<int> > I2);

void test_ari();

template<class T>
T sum(vector<T> & vec)
{
    T m = 0;
    for(int i=0; i<vec.size(); i++)
    {
        m = m + vec[i];
    }
    return m;
}

template<class T>
T mean(vector<T> & vec)
{
    T m = 0;
    for(int i=0; i<vec.size(); i++)
    {
        m = m + vec[i];
    }
    m = m / vec.size();
    return m;
}

template<class T>
T var(vector<T> & vec)
{
    T m = mean(vec);
    T var = 0;
    for(int i=0; i<vec.size(); i++)
    {
        var = var + (vec[i]-m)^2/vec.size();
    }
    return var;
}

}

namespace fileutils
{
vector<string> getdir (string dir);

bool file_exists (const std::string& name);
}
#endif // UTILS_H_INCLUDED
