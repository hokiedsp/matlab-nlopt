#include "mexObjectiveFunction.h"

#include <mexRuntimeError.h>

#include <algorithm>
#include <stdexcept>

#ifndef MX_CURRENT_API_VER
#define mxIsScalar(A) (mxGetNumberOfElements(A)==1)
#endif

mexObjectiveFunction::mexObjectiveFunction(nlopt_opt &optin, mxArray *mxFun, mxArray *mxHessMultFcn, mxArray *mxOutputFun)
    : prhs{mxFun, mxCreateDoubleMatrix(nlopt_get_dimension(optin), 1, mxREAL)}, // {fun, x}
      opt(optin),                                                               // nlopt_opt
      hessmult_args{NULL, NULL, NULL},     // {HessMultFcn x v}
      outfun_args{NULL, NULL, NULL, NULL}, // {OutputFun x optimValues state}
      lasterror(NULL), stop(false)
{
   // setup arguments for HessianMultFcn() evaluation only if function given
   if (!mxIsEmpty(mxHessMultFcn))
   {
      // initialize the arguments
      hessmult_args[0] = mxHessMultFcn;                                               // function handle
      hessmult_args[1] = mxCreateDoubleMatrix(nlopt_get_dimension(opt), 1, mxREAL); // x
      hessmult_args[2] = mxCreateDoubleMatrix(nlopt_get_dimension(opt), 1, mxREAL); // v
   }

   // setup arguments for OutputFun() evaluation only if function given
   if (!mxIsEmpty(mxOutputFun))
   {
      // initialize the arguments
      outfun_args[0] = mxOutputFun;                                // function handle
      outfun_args[1] = prhs[1];                                    // share the same mxArray with objective function
      outfun_args[2] = mexObjectiveFunction::create_optimvalues(); // optimvalues struct
   }
}

mexObjectiveFunction::~mexObjectiveFunction()
{
  mxDestroyArray(prhs[1]);
  outfun_args[1] = NULL; // same as prhs[1]
  for (int i = 1; i < 4; ++i)
    if (outfun_args[i])
      mxDestroyArray(outfun_args[i]);
  for (int i = 1; i < 3; ++i)
    if (hessmult_args[i])
      mxDestroyArray(hessmult_args[i]);
  if (lasterror)
    mxDestroyArray(lasterror);
}

double mexObjectiveFunction::obj_fun(unsigned n, const double *x, double *gradient, void *d_)
{
  mexObjectiveFunction &data = *(mexObjectiveFunction *)d_;
  double f = data.evalFun(n, x, gradient);
  if (data.stop)
    nlopt_force_stop(data.opt);
  return f;
};

void mexObjectiveFunction::precond_fun (unsigned n, const double *x, const double *v, double *vpre, void *f_data)
{
  mexObjectiveFunction &data = *(mexObjectiveFunction *)f_data;
  data.evalHessMultFcn(n, x, v, vpre);
  if (data.stop)
    nlopt_force_stop(data.opt);
};

double mexObjectiveFunction::evalFun(unsigned n, const double *x, double *gradient)
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

void mexObjectiveFunction::evalHessMultFcn(unsigned n, const double *x, const double *v, double *vpre)
{
   // prepare input and output arguments
   std::copy_n(x, n, mxGetPr(hessmult_args[1])); // copy the given x to input argument mxArray array
   std::copy_n(v, n, mxGetPr(hessmult_args[2])); // copy the given x to input argument mxArray array
   mxArray *plhs = NULL;

   // run the objective function
   if (call_matlab_feval_with_trap(1, &plhs, 3, hessmult_args) || // traps an error (incl. invalid # of arguments)
       !mxIsDouble(plhs) || mxIsComplex(plhs) || !(mxGetM(plhs) == 1 || mxGetN(plhs) == 1) || mxGetNumberOfElements(plhs) != n)
      stop = true;
   else
      std::copy_n(mxGetPr(plhs),n,vpre);

   if (plhs)
      mxDestroyArray(plhs);
}

mxArray *mexObjectiveFunction::create_optimvalues()
{
  const char *fields[22] = {"attainfactor", "cgiterations", "constrviolation", "degenerate", "directionalderivative", "firstorderopt",
                            "funccount", "fval", "gradient", "iteration", "lambda", "lssteplength", "maxfval", "positivedefinite", "procedure", "ratio",
                            "residual", "resnorm", "searchdirection", "stepaccept", "stepsize", "trustregionradius"};
  mxArray *rval = mxCreateStructMatrix(1, 1, 22, fields);
  for (int i = 0; i < 22; ++i)
    if (strcmp(fields[i], "funccount"))
      mxSetField(rval, 0, fields[i], mxCreateDoubleMatrix(0, 0, mxREAL));
    else
      mxSetField(rval, 0, fields[i], mxCreateDoubleScalar(0.0));

  return rval;
}

void mexObjectiveFunction::set_outfun_args(const char *state, mxArray *mxFval, mxArray *mxGradient)
{
  // update optimvalues.fval
  mxArray *field = mxGetField(outfun_args[2], 0, "fval");
  if (field)
    mxDestroyArray(field);
  mxSetField(outfun_args[2], 0, "fval", mxFval);

  // update optimvalues.gradient if gradient available
  field = mxGetField(outfun_args[2], 0, "gradient");
  if (mxGradient || !mxIsEmpty(field))
  {
    if (field)
      mxDestroyArray(field);
    mxSetField(outfun_args[2], 0, "gradient", mxGradient ? mxGradient : mxCreateDoubleMatrix(0, 0, mxREAL));
  }

  // update number of function evaluations
  *mxGetPr(mxGetField(outfun_args[2], 0, "funccount")) = nlopt_get_numevals(opt);

  // update the state (only do so if new state given)
  if (state)
  {
    if (outfun_args[3])
      mxDestroyArray(outfun_args[3]);
    outfun_args[3] = mxCreateString(state);
  }
}

mxArray *mexObjectiveFunction::evalGrad(mxArray *x)
{
  // assume output_args[1] already populated (last evalFun input)
  // set the optimvalues struct fields
  mxArray *plhs[2] = {NULL, NULL};
  prhs[1] = x;
  mxArray *MException = mexCallMATLABWithTrap(2, plhs, 2, prhs, "feval");
  if (MException) // maybe it does not compute gradient
  {
    mxDestroyArray(MException);
    return mxCreateDoubleMatrix(0, 0, mxREAL);
  }
  mxDestroyArray(plhs[0]);
  return plhs[1];
}

bool mexObjectiveFunction::evalOutputFun(bool init)
{
  // no OutputFun assigned, just return and continue
  if (!outfun_args[0]) return false;

  // assume output_args[1] already populated (last evalFun input)
  // set the optimvalues struct fields
  mxArray *plhs[2] = {NULL, NULL};
  mxArray *MException = mexCallMATLABWithTrap(2, plhs, 2, prhs, "feval");
  if (MException) // maybe it does not compute gradient
  {
    mxDestroyArray(MException);
    mexCallMATLAB(1, plhs, 2, prhs, "feval"); // try again just to retrieve the objective function value
  }

  // update optimvalues fields
  set_outfun_args(init ? "init" : "done", plhs[0], plhs[1]); // let go of the plhs content

  // run OutputFun
  mexCallMATLAB(1, plhs, 2, prhs, "feval");
  if (!((mxIsLogical(plhs[0])||mxIsNumeric(plhs[0])) && mxIsScalar(plhs[0])))
    throw mexRuntimeError("failedOutputFun","OutputFun must return a scalar logical.");

  // if init, change state to "iter" to prep for the forthcoming iterations
  if (init)
  {
    mxDestroyArray(outfun_args[3]);
    outfun_args[3] = mxCreateString("iter");
  }

  // return true if stop
  bool rval = (bool)mxGetScalar(plhs[0]);
  mxDestroyArray(plhs[0]);
  return rval;
}

bool mexObjectiveFunction::call_matlab_feval_with_trap(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
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

// void mexHessianFunction(unsigned n, const double *x, const double *v, double *vpre, void *d_)
// {
//   user_function_data *d = ((user_function_data *)d_)->dpre;
//   d->plhs[0] = d->plhs[1] = NULL;
//   memcpy(mxGetPr(d->prhs[d->xrhs]), x, n * sizeof(double));
//   memcpy(mxGetPr(d->prhs[d->xrhs + 1]), v, n * sizeof(double));

//   CHECK0(0 == mexCallMATLAB(1, d->plhs, d->nrhs, d->prhs, d->f),
//          "error calling user function");

//   CHECK0(mxIsDouble(d->plhs[0]) && !mxIsComplex(d->plhs[0]) && (mxGetM(d->plhs[0]) == 1 || mxGetN(d->plhs[0]) == 1) && mxGetM(d->plhs[0]) * mxGetN(d->plhs[0]) == n,
//          "vpre vector from user function is the wrong size");
//   memcpy(vpre, mxGetPr(d->plhs[0]), n * sizeof(double));
//   mxDestroyArray(d->plhs[0]);
//   d->neval++;
//   if (d->verbose)
//     mexPrintf("nlopt_optimize precond eval #%d\n", d->neval);
// }
