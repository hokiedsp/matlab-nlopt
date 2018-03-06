#pragma once

#include "mexObjectiveFunction.h"

#include <nlopt.h>
#include <mex.h>

#include <string>
#include <unordered_map>
#include <utility>

#define MEX_ACTION_ARGUMENTS const mxArray *mxObj, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The class that we are interfacing to
class mexNLopt
{
public:
  mexNLopt(const mxArray *mxObj, int nrhs, const mxArray *prhs[]) : opt(NULL) { init(prhs[0], prhs[1]); }

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
  ~mexNLopt() { if (opt) nlopt_destroy(opt); }

  /**
   * \brief Copy assignment
   */
  mexNLopt &operator=(const mexNLopt &src);

  /**
   * \brief Copy assignment
   */
  mexNLopt &operator=(mexNLopt &&src) noexcept;

  /**
   * \brief Returns the wrapping MATLAB class name
   */
  static std::string get_classname() { return "nlopt.options"; }; // must match the Matlab classname

  /**
   * \brief Performs static (object-independent) actions
   */
  static bool static_handler(std::string command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

  /**
   * \brief Performs object-dependent actions
   */
  bool action_handler(const mxArray *mxObj, const std::string &command, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

private:
  /**
 * \brief Return NLopt version string
 */
  static mxArray *getNLoptVersion();

  /**
 * \brief Return algorithm name given algorithm enum
 */
  static std::string get_algorithm_name_string(nlopt_algorithm a);

  /**
 * \brief Return algorithm name given algorithm enum
 */
  static mxArray *get_algorithm_name(nlopt_algorithm a);

  /**
 * \brief Return algorithm description given algorithm enum
 */
  static mxArray *get_algorithm_desc(nlopt_algorithm a);

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
  static mxArray *getAlgorithms(const mxArray *incl_desc);

  /**
 * \brief Returns algorithm given by mxString
 */
  static nlopt_algorithm find_algorithm_by_name(const mxArray *mxStr);

  nlopt_opt opt; /** nlopt_opt "object" (an opaque pointer) */

  /**
 * \brief Initialize nlopt object
 */
  void init(const mxArray *mxAlgorithm, const mxArray *mxDim);

  // GET ROUTINES

  void getAlgorithm(MEX_ACTION_ARGUMENTS) const { plhs[0] = mexNLopt::get_algorithm_name(nlopt_get_algorithm(opt)); }
  void getDimension(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar((double)nlopt_get_dimension(opt)); }

  /* stopping criteria: */
  void getStopVal(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_stopval(opt)); }
  void getFTolRel(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_ftol_rel(opt)); }
  void getXTolRel(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_xtol_rel(opt)); }
  void getFTolAbs(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_ftol_abs(opt)); }
  void getXTolAbs(MEX_ACTION_ARGUMENTS) const;
  void getMaxEval(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar((double)nlopt_get_maxeval(opt)); }
  void getMaxTime(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar(nlopt_get_maxtime(opt)); }
  void getPopulation(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar((double)nlopt_get_population(opt)); }
  void getVectorStorage(MEX_ACTION_ARGUMENTS) const { plhs[0] = mxCreateDoubleScalar((double)nlopt_get_vector_storage(opt)); }
  void getInitialStep(MEX_ACTION_ARGUMENTS) const;
  mxArray *getInitialStep(const mxArray *x0) const;

  // SET ROUTINES
  void setStopVal(MEX_ACTION_ARGUMENTS)
  {
    double sval = mxGetScalar(prhs[0]);
    nlopt_set_stopval(opt, isinf(sval) ? (sval > 0) ? HUGE_VAL : -HUGE_VAL : sval);
  }
  void setFTolRel(MEX_ACTION_ARGUMENTS) { nlopt_set_ftol_rel(opt, mxGetScalar(prhs[0])); }
  void setXTolRel(MEX_ACTION_ARGUMENTS) { nlopt_set_xtol_rel(opt, mxGetScalar(prhs[0])); }
  void setFTolAbs(MEX_ACTION_ARGUMENTS) { nlopt_set_ftol_abs(opt, mxGetScalar(prhs[0])); }
  void setXTolAbs(MEX_ACTION_ARGUMENTS)
  {
    if (mxIsScalar(prhs[0]))
      nlopt_set_xtol_abs1(opt, mxGetScalar(prhs[0]));
    else
      nlopt_set_xtol_abs(opt, mxGetPr(prhs[0]));
  }
  void setMaxEval(MEX_ACTION_ARGUMENTS) { nlopt_set_maxeval(opt, (int)mxGetScalar(prhs[0])); }
  void setMaxTime(MEX_ACTION_ARGUMENTS) { nlopt_set_maxtime(opt, mxGetScalar(prhs[0])); }
  void setPopulation(MEX_ACTION_ARGUMENTS) { nlopt_set_population(opt, (unsigned)mxGetScalar(prhs[0])); }
  void setVectorStorage(MEX_ACTION_ARGUMENTS) { nlopt_set_vector_storage(opt, (unsigned)mxGetScalar(prhs[0])); }
  void setInitialStep(MEX_ACTION_ARGUMENTS)
  {
    if (mxIsChar(prhs[0])) // 'auto'
      nlopt_set_initial_step(opt, NULL);
    else if (mxIsScalar(prhs[0]))
      nlopt_set_initial_step1(opt, mxGetScalar(prhs[0]));
    else
      nlopt_set_initial_step(opt, mxGetPr(prhs[0]));
  }

  // OPTIMIZATION ROUTINES
  void fminunc(MEX_ACTION_ARGUMENTS);
  void fmincon(MEX_ACTION_ARGUMENTS);
  void fminbnd(MEX_ACTION_ARGUMENTS);

  // common routine to run and gather the outcome
  static std::pair<nlopt_func,nlopt_precond> config_obj_fun(nlopt_opt opt, mexObjectiveFunction &data);
  static void run_n_report(int nlhs, mxArray *plhs[], nlopt_opt opt, mexObjectiveFunction &data);
  
  // routine to set subproblem options
  static void set_local_optimizer(nlopt_opt opt, const mxArray *mxObj);

  static void set_bounds(nlopt_opt opt, const mxArray *mxLB, const mxArray *mxUB);

  // map of actions
  typedef void (mexNLopt::*action_fcn)(MEX_ACTION_ARGUMENTS);
  static const std::unordered_map<std::string, action_fcn> action_map;

  typedef void (mexNLopt::*const_action_fcn)(MEX_ACTION_ARGUMENTS) const;
  static const std::unordered_map<std::string, const_action_fcn> const_action_map;
};
