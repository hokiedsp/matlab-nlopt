#pragma once

#include "mexObjectiveFunction.h"

#include <nlopt.h>
#include <mex.h>

struct mexConstraintFunction
{
  static double con_fun(unsigned n, const double *x, double *gradient, void *d_);
  static void vcon_fun(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data);

  bool isvector;
  nlopt_func fun;   // (unsigned n, const double *x, double *fc, void *func_data)
  nlopt_mfunc mfun; // (unsigned m, double *result, unsigned n, const double *x, double *fc, void *func_data)

  mxArray *prhs[2];          // feval mexMatlabCall input arguments for objective function evaluation
  
  nlopt_opt &opt;
  mxArray *&lasterror;        // mexObjectiveFunction's to store trapped MException 
  bool &stop;                 // mexObjectiveFunction's to flag if user issued a stop
  
  mexConstraintFunction(mxArray * mxFun, mexObjectiveFunction &data);
  ~mexConstraintFunction();
private:
  double evalFun(unsigned n, const double *x, double *gradient);
  void evalVecFun(unsigned m, unsigned n, const double* x, double *c, double* Cgrad);
  bool call_matlab_feval_with_trap(int nlhs, mxArray *plhs[], const int n, const double *x);
};
