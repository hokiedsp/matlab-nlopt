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
