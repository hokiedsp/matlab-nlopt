#include "mexConstraintFunction.h"

#include <mexRuntimeError.h>

#include <algorithm>
#include <stdexcept>

#ifndef MX_CURRENT_API_VER
#define mxIsScalar(A) (mxGetNumberOfElements(A)==1)
#endif

mexConstraintFunction::mexConstraintFunction(mxArray *mxFun, mexObjectiveFunction &data)
    : prhs{mxFun, mxCreateDoubleMatrix(nlopt_get_dimension(data.opt), 1, mxREAL)}, // {fun, x}
      opt(data.opt), lasterror(data.lasterror), stop(data.stop)
{
} 

mexConstraintFunction::~mexConstraintFunction()
{
  mxDestroyArray(prhs[1]);
}

double mexConstraintFunction::con_fun(unsigned n, const double *x, double *gradient, void *d_)
{
  mexConstraintFunction &data = *(mexConstraintFunction *)d_;
  double f = data.evalFun(n, x, gradient);
  if (data.stop)
    nlopt_force_stop(data.opt);
  return f;
};

void mexConstraintFunction::vcon_fun(unsigned m, double *result, unsigned n, const double *x, double *gradient, void *d_)
{
  mexConstraintFunction &data = *(mexConstraintFunction *)d_;
  data.evalVecFun(m, n, x, result, gradient);
  if (data.stop)
    nlopt_force_stop(data.opt);
};

double mexConstraintFunction::evalFun(unsigned n, const double *x, double *gradient)
{
  mxArray *plhs[2] = {NULL, NULL};

  // run the objective function
  if (call_matlab_feval_with_trap(gradient ? 2 : 1, plhs, n, x) ||       // trapped an error (incl. invalid # of arguments)
      !mxIsDouble(plhs[0]) || mxIsComplex(plhs[0]) || !mxIsScalar(plhs[0])) // must return a real double scalar value
    goto objfcn_failed;

  // prepare the function value to return
  double fval = mxGetScalar(plhs[0]);

  // assign the gradient vector to return if requested
  if (gradient)
  {
    // check for its validity first
    if (!mxIsDouble(plhs[1]) || mxIsComplex(plhs[1]) || !(mxGetM(plhs[1]) == 1 || mxGetN(plhs[1]) == 1) || mxGetNumberOfElements(plhs[1]) != n)
      goto objfcn_failed;

    std::copy_n(mxGetPr(plhs[1]), n, gradient);
  }

  mxDestroyArray(plhs[0]);
  mxDestroyArray(plhs[1]);

  return fval;

objfcn_failed:
   stop = true; // something went wrong, stop now
   if (plhs[0])
      mxDestroyArray(plhs[0]);
   if (plhs[1])
      mxDestroyArray(plhs[1]);
   return NAN;
}

void mexConstraintFunction::evalVecFun(unsigned m, unsigned n, const double* x, double *c, double* gradient)
{
  mxArray *plhs[2] = {NULL, NULL};

  // run the objective function
  if (call_matlab_feval_with_trap(gradient ? 2 : 1, plhs, n, x) ||       // trapped an error (incl. invalid # of arguments)
      !(mxGetM(plhs[0]) == 1 || mxGetN(plhs[0]) == 1) || mxGetNumberOfElements(plhs[0]) != n) // must return a real double scalar value
    goto objfcn_failed;

  // prepare the function value to return
  std::copy_n(mxGetPr(plhs[0]), n, c);

  // assign the gradient vector to return if requested
  if (gradient)
  {
    // check for its validity first
    if (!mxIsDouble(plhs[1]) || mxIsComplex(plhs[1]) || mxGetM(plhs[1]) != n || mxGetN(plhs[1]) != m)
      goto objfcn_failed;

    std::copy_n(mxGetPr(plhs[1]), n*m, gradient);
  }

  mxDestroyArray(plhs[0]);
  mxDestroyArray(plhs[1]);

  return;

objfcn_failed:
   stop = true; // something went wrong, stop now
   if (plhs[0])
      mxDestroyArray(plhs[0]);
   if (plhs[1])
      mxDestroyArray(plhs[1]);
}

bool mexConstraintFunction::call_matlab_feval_with_trap(int nlhs, mxArray *plhs[], const int n, const double *x)
{
    // prepare input and output arguments
  std::copy_n(x, n, mxGetPr(prhs[1])); // copy the given x to input argument mxArray array

  // call Matlab
  mxArray *new_error = mexCallMATLABWithTrap(nlhs, plhs, 2, prhs, "feval");

  // if errored out, keep the error
  if (new_error)
  {
    if (lasterror)
      mxDestroyArray(lasterror);
    lasterror = new_error;
    return true;
  }

  return false;
}
