function varargout = fminunc(fun,x0,options)
%NLOPT.FMINUNC   Find minimum of unconstrained multivariable function
% x = fminunc(fun,x0) starts at the point x0 and attempts to find a minimum
% x of the function described in fun. The point x0 can be a scalar, vector,
% or matrix.
%
%     Note:   Passing Extra Parameters explains how to pass extra
%     parameters to the objective function and nonlinear constraint
%     functions, if necessary.
%
%     fminunc is for nonlinear problems without constraints. If your
%     problem has constraints, generally use fmincon. See Optimization
%     Decision Table.
%
% x = fminunc(fun,x0,options) minimizes fun with the optimization options
% specified in options. Use optimoptions to set these options.
%
% [x,fval] = fminunc(...), for any syntax, returns the value of the
% objective function fun at the solution x.
%
% [x,fval,exitflag,output] = fminunc(...) additionally returns a value
% exitflag that describes the exit condition of fminunc, and a structure
% output with information about the optimization process.
%
% [x,fval,exitflag,output,grad,hessian] = fminunc(...) additionally
% returns:
%
%     grad — Gradient of fun at the solution x.
%
%     hessian — Hessian of fun at the solution x. See fminunc Hessian.

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
