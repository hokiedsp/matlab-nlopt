#pragma once

#include <nlopt.h>
#include <mex.h>

struct mexObjectiveFunction
{
  mxArray *prhs[2];          // feval mexMatlabCall input arguments for objective function evaluation
  nlopt_opt &opt;
  mxArray *hessmult_args[3];   // Hv = HessMultFcn(x,v)
  mxArray *outfun_args[4];   // feval mexMatlabCall input arguments for OutputFun function evaluation
  mxArray *lasterror;        // trapped MException
  bool stop;                 // true if user issued a stop
  
  mexObjectiveFunction(nlopt_opt & optin, mxArray * mxFun, mxArray * mxHessMultFcn, mxArray * mxOutputFun);
  ~mexObjectiveFunction();
  double evalFun(unsigned n, const double *x, double *gradient);
  void evalHessMultFcn(unsigned n, const double *x, const double *v, double *vpre);
  bool evalOutputFun(bool init);
  mxArray *evalGrad(mxArray *x);
private:
  static mxArray *create_optimvalues();
  void set_outfun_args(const char *state, mxArray *mxFval, mxArray *mxGradient = NULL);
  bool call_matlab_feval_with_trap(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);
};
