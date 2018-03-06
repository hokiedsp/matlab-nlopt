#pragma once

#include <nlopt.h>
#include <mex.h>

struct mexConstraintFunction
{
  bool isvector;
  nlopt_func fun;   // (unsigned n, const double *x, double *fc, void *func_data)
  nlopt_mfunc mfun; // (unsigned m, double *result, unsigned n, const double *x, double *fc, void *func_data)

  mxArray *prhs[2];          // feval mexMatlabCall input arguments for objective function evaluation
  mxArray *lasterror;        // trapped MException
  bool stop;                 // true if user issued a stop
  
  mexConstraintFunction(mxArray * mxFun, const nlopt_opt opt, const bool isvec);
  ~mexConstraintFunction();
private:
  double evalFun(unsigned n, const double *x, double *gradient);
  bool call_matlab_feval_with_trap(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);
};
