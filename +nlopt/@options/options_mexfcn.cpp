
#include "nlopt_algorigthm_idstr.h"

#include <mexObjectHandler.h>
#include <mex.h>

#include <nlopt.h>
#include <cctype>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

class mexNLopt;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mexObjectHandler<mexNLopt>(nlhs, plhs, nrhs, prhs);
}

// prototypes for the functions to be evaluated
struct mexObjectiveFunctionData { mxArray *prhs[2], *f; nlopt_opt &opt; };
static double mexObjectiveFunction(unsigned n, const double *x, double *gradient, void *d_);
static void mexObjectiveHessianFunction(unsigned n, const double *x, const double *v, double *vpre, void *d_);

// The class that we are interfacing to
class mexNLopt
{
public:
  mexNLopt(const mxArray *mxObj, int nrhs, const mxArray *prhs[]) : opt(NULL)
  {
    init(prhs[0], prhs[1]);
  }

  /**
   * \brief Copy constructor
   */
  mexNLopt(const mexNLopt &src) : opt(nlopt_copy(src.opt)) {}

  /**
   * \brief Move constructor
   */
  mexNLopt(mexNLopt &&src) : opt(std::move(src.opt)) {}

  /**
   * \brief Destructor
   */
  ~mexNLopt()
  {
    if (opt)
      nlopt_destroy(opt);
  }

  /**
   * \brief Copy assignment
   */
  mexNLopt &operator=(const mexNLopt &src)
  {
    if (opt)
      nlopt_destroy(opt);
    opt = nlopt_copy(src.opt);
    return *this;
  }

  /**
   * \brief Copy assignment
   */
  mexNLopt &operator=(mexNLopt &&src) noexcept
  {
    if (opt)
      nlopt_destroy(opt);
    opt = std::move(src.opt);
    return *this;
  }

  /**
   * \brief Returns the wrapping MATLAB class name
   */
  static std::string get_classname() { return "nlopt.options"; }; // must match the Matlab classname

  /**
   * \brief Performs static (object-independent) actions
   */
  static bool static_handler(std::string command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
  {
    if (command == "getNLoptVersion")
      plhs[0] = mexNLopt::getNLoptVersion();
    else if (command == "getAlgorithms")
      plhs[0] = mexNLopt::getAlgorithms(prhs[0]); // name_desc = getAlgorithms(incl_desc)
    else
      return false;
    return true;
  }

  /**
   * \brief Performs object-dependent actions
   */
  virtual bool action_handler(const mxArray *mxObj, const std::string &command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
  {
    // action command is looked up in the maps to get the pointer to the member function
    //    action_map - unqualified member functions
    //    const_action_map - const-qualified member functions
    // maps are static const variable, defined at the bottom of this file
    auto action_iterator = action_map.find(command);
    if (action_iterator != action_map.end())
      (this->*(action_iterator->second))(nlhs, plhs, nrhs, prhs);
    else
    {
      auto const_action_iterator = const_action_map.find(command);
      if (const_action_iterator == const_action_map.end())
        return false;
      (this->*(const_action_iterator->second))(nlhs, plhs, nrhs, prhs);
    }
    return true;
  }

private:
  /**
 * \brief Return NLopt version string
 */
  static mxArray *getNLoptVersion()
  {
    int major, minor, bugfix;
    nlopt_version(&major, &minor, &bugfix);

    std::ostringstream buf;
    buf << major << '.' << minor << '.' << bugfix;

    return mxCreateString(buf.str().c_str());
  }

  /**
 * \brief Return algorithm name given algorithm enum
 */
  static std::string get_algorithm_name_string(nlopt_algorithm a)
  {
    if (a == NLOPT_NUM_ALGORITHMS)
      throw mexRuntimeError(mexNLopt::get_classname() + ":invalidAlgorithm", "Invalid algorithm name.");
    std::string name = nlopt_algorithm_idstrs[(int)a] + 6;
    std::for_each(name.begin(), name.end(), [](auto &ch) { ch = tolower(ch); });
    return name;
  }
  
  /**
 * \brief Return algorithm name given algorithm enum
 */
  static mxArray *get_algorithm_name(nlopt_algorithm a)
  {
    return mxCreateString(get_algorithm_name_string(a).c_str());
  }

  /**
 * \brief Return algorithm description given algorithm enum
 */
  static mxArray *get_algorithm_desc(nlopt_algorithm a)
  {
    if (a == NLOPT_NUM_ALGORITHMS)
      throw mexRuntimeError(mexNLopt::get_classname() + ":invalidAlgorithm", "Invalid algorithm name.");
    return mxCreateString(nlopt_algorithm_name(a));
  }

  /**
 * \brief Iterate over algorithms performing unary function on each element
 * 
 * \note signature for f: \code bool f(const nlopt_algorithm a) \endcode 
 * \note if f returns false, iteration stops immediately
 */
  template <class UnaryFunction>
  static void for_each_algorithm(UnaryFunction f)
  {
    for (int n = 0; n < (int)NLOPT_NUM_ALGORITHMS && f((nlopt_algorithm)n); ++n)
      ;
  }

  /**
 * \brief Returns a cellstr array of all available algorithms
 */
  static mxArray *getAlgorithms(const mxArray *incl_desc)
  {
    int N = (int)NLOPT_NUM_ALGORITHMS;
    mxArray *rval = mxCreateCellMatrix(N, incl_desc ? 2 : 1);
    for_each_algorithm([&](const nlopt_algorithm a) -> bool {
      mxSetCell(rval, (mwIndex)a, mexNLopt::get_algorithm_name(a));
      if (mxGetLogicals(incl_desc)[0])
        mxSetCell(rval, N + (mwIndex)a, mexNLopt::get_algorithm_desc(a));
      return true;
    });
    return rval;
  }

  /**
 * \brief Returns algorithm given by mxString
 */
  static nlopt_algorithm find_algorithm_by_name(const mxArray *mxStr)
  {
    std::string name = mexGetString(mxStr);
    nlopt_algorithm algo = NLOPT_NUM_ALGORITHMS;
    for_each_algorithm([&](const nlopt_algorithm a) -> bool {
      bool rval;
      if (rval = (name == get_algorithm_name_string(a)))
        algo = a;
      return !rval;
    });

    // not found
    if (algo == NLOPT_NUM_ALGORITHMS)
      throw mexRuntimeError(mexNLopt::get_classname() + ":invalidAlgorithm", "Invalid algorithm name.");

    return algo;
  }

  // typedef double (*nlopt_func)(unsigned n, const double *x,
  // 			     double *gradient, /* NULL if not needed */
  // 			     void *func_data);

  // typedef void (*nlopt_mfunc)(unsigned m, double *result,
  // 			    unsigned n, const double *x,
  // 			     double *gradient, /* NULL if not needed */
  // 			     void *func_data);

  // /* A preconditioner, which preconditions v at x to return vpre.
  //    (The meaning of "preconditioning" is algorithm-dependent.) */
  // typedef void (*nlopt_precond)(unsigned n, const double *x, const double *v,
  // 			      double *vpre, void *data);

  // typedef enum {
  //      NLOPT_FAILURE = -1, /* generic failure code */
  //      NLOPT_INVALID_ARGS = -2,
  //      NLOPT_OUT_OF_MEMORY = -3,
  //      NLOPT_ROUNDOFF_LIMITED = -4,
  //      NLOPT_FORCED_STOP = -5,
  //      NLOPT_SUCCESS = 1, /* generic success code */
  //      NLOPT_STOPVAL_REACHED = 2,
  //      NLOPT_FTOL_REACHED = 3,
  //      NLOPT_XTOL_REACHED = 4,
  //      NLOPT_MAXEVAL_REACHED = 5,
  //      NLOPT_MAXTIME_REACHED = 6
  // } nlopt_result;

  nlopt_opt opt; /** nlopt_opt "object" (an opaque pointer) */

  /**
 * \brief Initialize nlopt object
 */
  void init(const mxArray *mxAlgorithm, const mxArray *mxDim)
  {
    nlopt_algorithm alg = find_algorithm_by_name(mxAlgorithm); // throw exception if not a valid algorithm name
    int dim = (int)mxGetScalar(mxDim);                         // assume prevalidated in MATLAB

    if (opt)
      nlopt_destroy(opt);

    opt = nlopt_create(alg, dim);
    if (opt == NULL)
      throw mexRuntimeError(mexNLopt::get_classname() + ":failedNLoptCreation", "nlopt_create() failed to create a new nlopt object.");
  }

  // GET ROUTINES

  void getAlgorithm(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mexNLopt::get_algorithm_name(nlopt_get_algorithm(opt)); }
  void getDimension(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar((double)nlopt_get_dimension(opt)); }

  /* stopping criteria: */
  void getStopVal(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_stopval(opt)); }
  void getFTolRel(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_ftol_rel(opt)); }
  void getXTolRel(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_xtol_rel(opt)); }
  void getFTolAbs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_ftol_abs(opt)); }
  void getXTolAbs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const
  {
    unsigned d = nlopt_get_dimension(opt);
    std::vector<double> tol(d);
    if (nlopt_get_xtol_abs(opt, tol.data()) < 0)
      throw mexRuntimeError(mexNLopt::get_classname() + ":failedNLoptGet", "nlopt_get_xtol_abs() failed.");
    double tol0 = tol[0];
    if (d == 1 || std::all_of(tol.begin() + 1, tol.end(), [tol0](const double t) { return t == tol0; }))
      plhs[0] = mxCreateDoubleScalar(tol0);
    else
    {
      mxArray *rval = mxCreateDoubleMatrix(d, 1, mxREAL);
      std::copy(tol.begin(), tol.end(), mxGetPr(rval));
      plhs[0] = rval;
    }
  }
  void getMaxEval(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar((double)nlopt_get_maxeval(opt)); }
  void getMaxTime(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_maxtime(opt)); }
  void getPopulation(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar((double)nlopt_get_population(opt)); }
  void getVectorStorage(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const { plhs[0] = mxCreateDoubleScalar((double)nlopt_get_vector_storage(opt)); }
  void getInitialStep(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const
  {
    unsigned d = nlopt_get_dimension(opt);
    std::vector<double> steps(d);
    std::vector<double> x0(d, 1.0);
    if (nlopt_get_initial_step(opt, x0.data(), steps.data()) < 0)
      throw mexRuntimeError(mexNLopt::get_classname() + ":failedNLoptGet", "nlopt_get_initial_step() failed.");
    double dx0 = steps[0];
    if (d == 1 || std::all_of(steps.begin() + 1, steps.end(), [dx0](const double dx) { return dx == dx0; }))
      plhs[0] = mxCreateDoubleScalar(dx0);
    else
    {
      mxArray *rval = mxCreateDoubleMatrix(d, 1, mxREAL);
      std::copy(steps.begin(), steps.end(), mxGetPr(rval));
      plhs[0] = rval;
    }
  }

  mxArray *getInitialStep(const mxArray *x0) const
  {
    unsigned d = nlopt_get_dimension(opt);
    std::vector<double> steps(d);
    if (nlopt_get_initial_step(opt, mxGetPr(x0), steps.data()) < 0)
      throw mexRuntimeError(mexNLopt::get_classname() + ":failedNLoptGet", "nlopt_get_initial_step() failed.");
    double dx0 = steps[0];
    if (d == 1 || std::all_of(steps.begin() + 1, steps.end(), [dx0](const double dx) { return dx == dx0; }))
      return mxCreateDoubleScalar(dx0);
    else
    {
      mxArray *rval = mxCreateDoubleMatrix(d, 1, mxREAL);
      std::copy(steps.begin(), steps.end(), mxGetPr(rval));
      return rval;
    }
  }

  // SET ROUTINES
  void setStopVal(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
  {
    double sval = mxGetScalar(prhs[0]);
    nlopt_set_stopval(opt, isinf(sval) ? (sval > 0) ? HUGE_VAL : -HUGE_VAL : sval);
  }
  void setFTolRel(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { nlopt_set_ftol_rel(opt, mxGetScalar(prhs[0])); }
  void setXTolRel(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { nlopt_set_xtol_rel(opt, mxGetScalar(prhs[0])); }
  void setFTolAbs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { nlopt_set_ftol_abs(opt, mxGetScalar(prhs[0])); }
  void setXTolAbs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
  {
    if (mxIsScalar(prhs[0]))
      nlopt_set_xtol_abs1(opt, mxGetScalar(prhs[0]));
    else
      nlopt_set_xtol_abs(opt, mxGetPr(prhs[0]));
  }
  void setMaxEval(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { nlopt_set_maxeval(opt, (int)mxGetScalar(prhs[0])); }
  void setMaxTime(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { nlopt_set_maxtime(opt, mxGetScalar(prhs[0])); }
  void setPopulation(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { nlopt_set_population(opt, (unsigned)mxGetScalar(prhs[0])); }
  void setVectorStorage(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { nlopt_set_vector_storage(opt, (unsigned)mxGetScalar(prhs[0])); }
  void setInitialStep(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
  {
    if (mxIsChar(prhs[0])) // 'auto'
      nlopt_set_initial_step(opt, NULL);
    else if (mxIsScalar(prhs[0]))
      nlopt_set_initial_step1(opt, mxGetScalar(prhs[0]));
    else
      nlopt_set_initial_step(opt, mxGetPr(prhs[0]));
  }

  // OPTIMIZATION ROUTINES
  void fminunc(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
  { // x = fminunc(fun,x0)
    unsigned d = nlopt_get_dimension(opt);
    mexObjectiveFunctionData data = {{(mxArray*) prhs[0], mxCreateDoubleMatrix(d,1,mxREAL)}, // fun, x
                                     NULL, // f
                                     opt}; // nlopt_opt
    const mxArray *mxFun = prhs[0]; // prevalidated function_handle
    if (false) // if Hessian given
    {
      // CHECK(mxIsChar(mx) || mxIsFunctionHandle(mx),
      //       "pre must contain function handles or function names");
      //   dpre.prhs[0] = mx;
      //   strcpy(dpre.f, "feval");
      //   dpre.nrhs = 3;
      //   dpre.xrhs = 1;
      // dpre.verbose = d.verbose > 2;
      // dpre.opt = opt;
      // dpre.neval = 0;
      // dpre.prhs[dpre.xrhs] = d.prhs[d.xrhs];
      // dpre.prhs[d.xrhs + 1] = mxCreateDoubleMatrix(1, n, mxREAL);
      // d.dpre = &dpre;
      //   nlopt_set_precond_min_objective(opt, user_function, user_pre, &d);
    }
    else // else just objective function
    {
      nlopt_set_min_objective(opt, mexObjectiveFunction, &data);
    }

    plhs[0] = mxDuplicateArray(prhs[1]); // create output x vector from prevalidated initial x
    double *x = mxGetPr(plhs[0]);   // get the vector
    double fval;
    if (nlopt_optimize(opt, x, &fval)<0)
    {
      throw mexRuntimeError(mexNLopt::get_classname() + ":fminuncFailed", "nlopt_optimize() failed.");
    }
  }

  // NLOPT_EXTERN(nlopt_result) nlopt_set_local_optimizer(nlopt_opt opt,
  // 						    const nlopt_opt local_opt);

  // NLOPT_EXTERN(nlopt_result) nlopt_set_default_initial_step(nlopt_opt opt,
  // 							 const double *x);

  // algorithm outcome
  // NLOPT_EXTERN(int) nlopt_get_numevals(const nlopt_opt opt);

  // NLOPT_EXTERN(nlopt_result) nlopt_force_stop(nlopt_opt opt);
  // NLOPT_EXTERN(nlopt_result) nlopt_set_force_stop(nlopt_opt opt, int val);
  // NLOPT_EXTERN(int) nlopt_get_force_stop(const nlopt_opt opt);

  /* more algorithm-specific parameters */

  /* the only immutable parameters of an optimization are the algorithm and
   the dimension n of the problem, since changing either of these could
   have side-effects on lots of other parameters */

  // NLOPT_EXTERN(nlopt_result) nlopt_optimize(nlopt_opt opt, double *x,
  // 					 double *opt_f);

  // NLOPT_EXTERN(nlopt_result) nlopt_set_min_objective(nlopt_opt opt, nlopt_func f,
  // 						  void *f_data);
  // NLOPT_EXTERN(nlopt_result) nlopt_set_precond_min_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);

  // NLOPT_EXTERN(nlopt_result) nlopt_set_max_objective(nlopt_opt opt, nlopt_func f,
  // 						  void *f_data);
  // NLOPT_EXTERN(nlopt_result) nlopt_set_precond_max_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre, void *f_data);

  // NLOPT_EXTERN(const char*) nlopt_get_errmsg(nlopt_opt opt);

  /* constraints: */

  // NLOPT_EXTERN(nlopt_result) nlopt_set_lower_bounds(nlopt_opt opt,
  // 						 const double *lb);
  // NLOPT_EXTERN(nlopt_result) nlopt_set_lower_bounds1(nlopt_opt opt, double lb);
  // NLOPT_EXTERN(nlopt_result) nlopt_get_lower_bounds(const nlopt_opt opt,
  // 						 double *lb);
  // NLOPT_EXTERN(nlopt_result) nlopt_set_upper_bounds(nlopt_opt opt,
  // 						 const double *ub);
  // NLOPT_EXTERN(nlopt_result) nlopt_set_upper_bounds1(nlopt_opt opt, double ub);
  // NLOPT_EXTERN(nlopt_result) nlopt_get_upper_bounds(const nlopt_opt opt,
  // 						 double *ub);

  // NLOPT_EXTERN(nlopt_result) nlopt_remove_inequality_constraints(nlopt_opt opt);
  // NLOPT_EXTERN(nlopt_result) nlopt_add_inequality_constraint(nlopt_opt opt,
  // 							  nlopt_func fc,
  // 							  void *fc_data,
  // 							  double tol);
  // NLOPT_EXTERN(nlopt_result) nlopt_add_precond_inequality_constraint(
  //      nlopt_opt opt, nlopt_func fc, nlopt_precond pre, void *fc_data,
  //      double tol);
  // NLOPT_EXTERN(nlopt_result) nlopt_add_inequality_mconstraint(nlopt_opt opt,
  // 							    unsigned m,
  // 							    nlopt_mfunc fc,
  // 							    void *fc_data,
  // 							    const double *tol);

  // NLOPT_EXTERN(nlopt_result) nlopt_remove_equality_constraints(nlopt_opt opt);
  // NLOPT_EXTERN(nlopt_result) nlopt_add_equality_constraint(nlopt_opt opt,
  // 							nlopt_func h,
  // 							void *h_data,
  // 							double tol);
  // NLOPT_EXTERN(nlopt_result) nlopt_add_precond_equality_constraint(
  //      nlopt_opt opt, nlopt_func h, nlopt_precond pre, void *h_data,
  //      double tol);
  // NLOPT_EXTERN(nlopt_result) nlopt_add_equality_mconstraint(nlopt_opt opt,
  // 							  unsigned m,
  // 							  nlopt_mfunc h,
  // 							  void *h_data,
  // 							  const double *tol);

  // map of actions
  typedef void (mexNLopt::*action_fcn)(int, mxArray *[], int, const mxArray *[]);
  static const std::unordered_map<std::string, action_fcn> action_map;

  typedef void (mexNLopt::*const_action_fcn)(int, mxArray *[], int, const mxArray *[]) const;
  static const std::unordered_map<std::string, const_action_fcn> const_action_map;
};

const std::unordered_map<std::string, mexNLopt::action_fcn> mexNLopt::action_map =
    {
        {"fminunc", &mexNLopt::fminunc},
        {"setFunctionStopValue", &mexNLopt::setStopVal},
        {"setFunctionAbsoluteTolerance", &mexNLopt::setFTolAbs},
        {"setMaxFunctionEvaluations", &mexNLopt::setMaxEval},
        {"setFunctionRelativeTolerance", &mexNLopt::setFTolRel},
        {"setInitialStepSize", &mexNLopt::setInitialStep},
        {"setPopulation", &mexNLopt::setPopulation},
        {"setStepAbsoluteTolerance", &mexNLopt::setXTolAbs},
        {"setStepRelativeTolerance", &mexNLopt::setXTolRel},
        {"setMaxEvaluationDuration", &mexNLopt::setMaxTime},
        {"setVectorStorage", &mexNLopt::setVectorStorage}};

const std::unordered_map<std::string, mexNLopt::const_action_fcn> mexNLopt::const_action_map =
    {
        {"getAlgorithm", &mexNLopt::getAlgorithm},
        {"getDimension", &mexNLopt::getDimension},
        {"getFunctionAbsoluteTolerance", &mexNLopt::getFTolAbs},
        {"getFunctionRelativeTolerance", &mexNLopt::getFTolRel},
        {"getFunctionStopValue", &mexNLopt::getStopVal},
        {"getInitialStepSize", &mexNLopt::getInitialStep},
        {"getMaxEvaluationDuration", &mexNLopt::getMaxTime},
        {"getMaxFunctionEvaluations", &mexNLopt::getMaxEval},
        {"getPopulation", &mexNLopt::getPopulation},
        {"getStepAbsoluteTolerance", &mexNLopt::getXTolAbs},
        {"getStepRelativeTolerance", &mexNLopt::getXTolRel},
        {"getVectorStorage", &mexNLopt::getVectorStorage}};

static double mexObjectiveFunction(unsigned n, const double *x, double *gradient, void *data_)
{
  mexObjectiveFunctionData &data = *(mexObjectiveFunctionData *)data_;
  mxArray *plhs[2] = {NULL, NULL};
  bool failed = false;

  // copy the given x to input argument mxArray array
  std::copy_n(x, n, mxGetPr(data.prhs[1]));
  
  mxArray *MException = mexCallMATLABWithTrap(gradient ? 2 : 1, plhs, 2, data.prhs, "feval");
  if (MException || !mxIsDouble(plhs[0]) || mxIsComplex(plhs[0]) || !mxIsScalar(plhs[0])) // trapped an error
  {
    failed = true;
    goto objfcn_end;
  }

  double f = mxGetScalar(plhs[0]);
  if (data.f) mxDestroyArray(data.f);
  data.f = plhs[0];

  if (gradient)
  {
    if (!mxIsDouble(plhs[1]) || mxIsComplex(plhs[1]) || !(mxGetM(plhs[1]) == 1 || mxGetN(plhs[1]) == 1) || mxGetNumberOfElements(plhs[1])!=n)
    {
      failed = true;
      goto objfcn_end;
    }
    
    std::copy_n(mxGetPr(plhs[1]), n, gradient);
    mxDestroyArray(plhs[1]);
  }

objfcn_end:
  // if (d->verbose)
  //   mexPrintf("nlopt_optimize eval #%d: %g\n", d->neval, f);
  if (failed || mxIsNaN(f))
    nlopt_force_stop(data.opt);
  return f;
}

// static void mexHessianFunction(unsigned n, const double *x, const double *v, double *vpre, void *d_)
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
