#include <fstream>

using namespace std;

void main()
{
    string outfile = "C:/Users/tikum/Documents/GitHub/matlab-nlopt/build/+nlopt/@options/nlopt_algorigthm_idstr.h";
    string infile = "C:/Users/tikum/Documents/GitHub/matlab-nlopt/build/nlopt/include/nlopt.h";

    ofstream out(outfile);
    out << "#pragma once" << endl
    <<"#include ""nlopt.h""" << endl
    << endl
    <<"static const char nlopt_algorithm_idstrs[NLOPT_NUM_ALGORITHMS][256] = {" << endl;

    std::ifstream in(infile);

}
// static const char nlopt_algorithm_short_names[NLOPT_NUM_ALGORITHMS][256] = {

// nlopt.h
// find:      "NLOPT_GN_DIRECT = 0,"
// then find "} nlopt_algorithm;"
// strip "= 0"
//
// detect any "#ifdef NLOPT_CXX"..."#else"..."#endif" and only use else
// wrap each name in a pair of double quotes
// 
//}
