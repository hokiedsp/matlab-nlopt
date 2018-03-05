#include "mexObjectiveFunctionData.h"

#include <algorithm>
#include <stdexcept>

double mexObjectiveFunction(unsigned n, const double *x, double *gradient, void *data_)
{
  mexObjectiveFunctionData &data = *(mexObjectiveFunctionData *)data_;
  double f = data.evalFun(n, x, gradient);

  if (data.stop)
    nlopt_force_stop(data.opt);
  
  // if (d->verbose)
  //   mexPrintf("nlopt_optimize eval #%d: %g\n", d->neval, f);
    
  return f;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

mexObjectiveFunctionData::mexObjectiveFunctionData(nlopt_opt &optin, mxArray *mxFun, mxArray *mxOutputFun)
    : prhs{mxFun, mxCreateDoubleMatrix(nlopt_get_dimension(optin), 1, mxREAL)}, // {fun, x}
      opt(optin),                                                                 // nlopt_opt
      outfun_args{NULL, NULL, NULL, NULL},                                      // {OutputFun x optimValues state}
      lasterror(NULL), stop(false)
{
  // if no output function, nothing to configure
  if (!mxIsEmpty(mxOutputFun))
  {
    // initialize the arguments
    outfun_args[0] = mxOutputFun;                                    // function handle
    outfun_args[1] = prhs[1];                                        // share the same mxArray with objective function
    outfun_args[2] = mexObjectiveFunctionData::create_optimvalues(); // optimvalues struct
  }
}

mexObjectiveFunctionData::~mexObjectiveFunctionData()
{
  mxDestroyArray(prhs[1]);
  outfun_args[1] = NULL; // same as prhs[1]
  for (int i = 1; i < 4; ++i)
    if (outfun_args[i])
      mxDestroyArray(outfun_args[i]);
  if (lasterror)
    mxDestroyArray(lasterror);
}

double mexObjectiveFunctionData::evalFun(unsigned n, const double *x, double *gradient)
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
    mexPrintf("x = [%f,%f], f = %f\n",mxGetPr(outfun_args[1])[0],mxGetPr(outfun_args[1])[1],mxGetScalar(plhs[0]));

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
  if (plhs[0])
    mxDestroyArray(plhs[0]);
  if (plhs[1])
    mxDestroyArray(plhs[1]);
  return NAN;
}

mxArray *mexObjectiveFunctionData::create_optimvalues()
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

void mexObjectiveFunctionData::set_outfun_args(const char *state, mxArray *mxFval, mxArray *mxGradient)
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

mxArray *mexObjectiveFunctionData::evalGrad(mxArray *x)
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

bool mexObjectiveFunctionData::evalOutputFun(bool init)
{
  // no OutputFun assigned, just return and continue
  if (!outfun_args[0]) return false;

  // // assume output_args[1] already populated (last evalFun input)
  // // set the optimvalues struct fields
  // mxArray *plhs[2] = {NULL, NULL};
  // mxArray *MException = mexCallMATLABWithTrap(2, plhs, 2, prhs, "feval");
  // if (MException) // maybe it does not compute gradient
  // {
  //   mxDestroyArray(MException);
  //   mexCallMATLAB(1, plhs, 2, prhs, "feval"); // try again just to retrieve the objective function value
  // }

  // // update optimvalues fields
  // set_outfun_args(init ? "init" : "done", plhs[0], plhs[1]); // let go of the plhs content

  // // run OutputFun
  // mexCallMATLAB(1, plhs, 2, prhs, "feval");
  // if (!(mxIsLogical(plhs[0]) && mxIsScalar(plhs[0])))
  //   throw std::runtime_error("OutputFun must return a scalar logical.");

  // // if init, change state to "iter" to prep for the forthcoming iterations
  // if (init)
  // {
  //   mxDestroyArray(outfun_args[3]);
  //   outfun_args[3] = mxCreateString("iter");
  // }

  // // return true if stop
  // bool rval = *mxGetLogicals(plhs[0]);
  // mxDestroyArray(plhs[0]);
  // return rval;
  return false;
}

bool mexObjectiveFunctionData::call_matlab_feval_with_trap(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
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
