#pragma once

#include <nlopt.h>
#include <mex.h>

double mexObjectiveFunction(unsigned n, const double *x, double *gradient, void *d_);
// void mexHessianFunction(unsigned n, const double *x, const double *v, double *vpre, void *d_);

struct mexObjectiveFunctionData
{
  mxArray *prhs[2];          // feval mexMatlabCall input arguments for objective function evaluation
  nlopt_opt &opt;
  mxArray *outfun_args[4];   // feval mexMatlabCall input arguments for OutputFun function evaluation
  mxArray *lasterror;        // trapped MException
  bool stop;                 // true if user issued a stop
  
  mexObjectiveFunctionData(nlopt_opt & optin, mxArray * mxFun, mxArray * mxOutputFun);
  ~mexObjectiveFunctionData();
  double evalFun(unsigned n, const double *x, double *gradient);
  bool evalOutputFun(bool init);
  mxArray *evalGrad(mxArray *x);
private:
  static mxArray *create_optimvalues();
  void set_outfun_args(const char *state, mxArray *mxFval, mxArray *mxGradient = NULL);
  bool call_matlab_feval_with_trap(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);
};
