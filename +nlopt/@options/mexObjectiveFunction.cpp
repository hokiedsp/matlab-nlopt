#include "mexConstraintFunction.h"

#include <mexRuntimeError.h>

#include <algorithm>
#include <stdexcept>

mexConstraintFunction::mexConstraintFunction(mxArray *mxFun, const nlopt_opt opt, const bool isvec)
    : isvector(isvec),
      prhs{mxFun, mxCreateDoubleMatrix(nlopt_get_dimension(optin), 1, mxREAL)}, // {fun, x}
      lasterror(NULL), stop(false)
{
  if (isvec)
    fun = [](unsigned n, const double *x, double *fc, void *d_) -> double {
      mexConstraintFunction &data = *(mexConstraintFunction *)d_;
      double f = data.evalFun(n, x, fc);
      if (data.stop)
        nlopt_force_stop(data.opt);
      return f;
    };
  else
    mfun = [](unsigned m, double *result, unsigned n, const double *x, double *fc, void *d_) -> double {
      mexConstraintFunction &data = *(mexConstraintFunction *)d_;
      double f = data.evalFun(n, x, gradient);
      if (data.stop)
        nlopt_force_stop(data.opt);
      return f;
    };
}

mexConstraintFunction::~mexConstraintFunction()
{
  mxDestroyArray(prhs[1]);
}

double mexConstraintFunction::evalFun(unsigned n, const double *x, double *gradient)
{
  // prepare input and output arguments
  std::copy_n(x, n, mxGetPr(prhs[1])); // copy the given x to input argument mxArray array
  mxArray *plhs[2] = {NULL, NULL};

  // run the objective function
  if (call_matlab_feval_with_trap(gradient ? 2 : 1, plhs, 2, prhs) ||       // trapped an error (incl. invalid # of arguments)
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

  // run OutputFun if specified
  if (outfun_args[0])
  {
    // update optimvalues
    set_outfun_args("", plhs[0], plhs[1]); // move the evalutaed values to outfun_args
    plhs[0] = plhs[1] = NULL;

    // call OutputFun
    if (call_matlab_feval_with_trap(1, plhs, 4, outfun_args) ||
        !(mxIsLogical(plhs[0]) && mxIsScalar(plhs[0])))
      goto objfcn_failed;

    // return true if stop
    stop = *mxGetLogicals(plhs[0]);
    mxDestroyArray(plhs[0]);
  }
  else
  {
    mxDestroyArray(plhs[0]);
    mxDestroyArray(plhs[1]);
  }

  return fval;

objfcn_failed:
   stop = true; // something went wrong, stop now
   if (plhs[0])
      mxDestroyArray(plhs[0]);
   if (plhs[1])
      mxDestroyArray(plhs[1]);
   return NAN;
}


bool mexConstraintFunction::call_matlab_feval_with_trap(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  // call Matlab
  mxArray *new_error = mexCallMATLABWithTrap(nlhs, plhs, nrhs, prhs, "feval");

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
