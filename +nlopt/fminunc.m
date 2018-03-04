function varargout = fminunc(fun,x0,options)
%NLOPT.FMINUNC   Find minimum of unconstrained multivariable function
%   x = nlopt.fminunc(fun,x0) starts at the point x0 and attempts to find a
%   minimum x of the objective function given as a function_handle in fun.
%   The point x0 can be a scalar or vector.
%
%   The objective function must have the form:
%
%      f = fun(x) or
%      [f,df] = fun(x)
%
%   where f is the function value evaluated at x and df is the gradient
%   of the function at x. If the gradient is provided, it uses 'ld_lbfgs'
%   algorithm. If the derivative is not provided, it uses 'ld_neldermead'
%   algorithm. See nlopt.getAlgorithms() to get the complete list of
%   available algorithms.
%
%   x = nlopt.fminunc(fun,x0,options) minimizes fun with the optimization
%   options specified in options. Use nlopt.options to set these options.
%
%   [x,fval] = fminunc(...) returns the value of the objective function fun
%   at the solution x.
%
%   [x,fval,exitflag,output] = fminunc(...) additionally returns a value
%   exitflag that describes the exit condition of fminunc, and a structure
%   output with information about the optimization process.
%
%   [x,fval,exitflag,output,grad,hessian] = fminunc(...) additionally
%   returns:
%
%      grad    - Gradient of fun at the solution x. Empty if FUN does not
%                return gradient
%      hessian - Hessian of fun at the solution x. Empty if HessianFun is
%                not given.
%
%   See also: nlopt.getAlgorithms, nlopt.options, nlopt.fmincon

narginchk(2,3);
nargoutchk(0,6);
if ~isa(fun,'function_handle')
   error('FUN must be a function handle.');
end
if nargin>2
   validateattributes(options,{'nlopt.options'},{'scalar'},'nlopt.fminunc','options');
   validateattributes(x0,{'double'},{'vector','numel',options.Dimension,'real','finite'});
else
   validateattributes(x0,{'double'},{'vector','real','finite'});
   try
      [~,~] = fun(x0);
      options = nlopt.options('ld_lbfgs',numel(x0));
   catch
      options = nlopt.options('ln_neldermead',numel(x0));
   end
end

[varargout{1:nargout}] = options.fminunc(fun,x0);
