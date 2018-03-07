function varargout = fminbnd(fun,lb,ub,options)
%nlopt.fminbnd Find bounded minimum of a multivariable function
%   nlopt.fminbnd(FUN,x1,x2) attempts to find a local minimizer X of the
%   function FUN in the interval x1 <= X <= x2.  FUN is a function handle.
%   FUN accepts input X and returns a scalar function value F
%   evaluated at X.
%
%   X = nlopt.fminbnd(FUN,x1,x2,OPTIONS) minimizes with the specified
%   algorithm and its options, specified by the nlopt.options object in
%   OPTIONS. Only the algorithms which supports bounded constraints would
%   work for this function and use a global optimization algorithm as the
%   initial values are not given. If an approximate location is given, use
%   nlopt.fmincon() function with a local optimization algorithm.
%
%   [X,FVAL] = nlopt.fminbnd(...) also returns the value of the objective
%   function, FVAL, computed in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = nlopt.fminbnd(...) also returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%
%    1  nlopt.fminbnd converged with a solution X based on OPTIONS.TolX. 0
%    Maximum number of function evaluations or iterations reached.
%   -1  Algorithm terminated by the output function. -2  Bounds are
%   inconsistent (that is, ax > bx).
%
%   [X,FVAL,EXITFLAG,OUTPUT] = nlopt.fminbnd(...) also returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   See also nlopt.options, nlopt.fminbnd, nlopt.fminunc

narginchk(3,4);
nargoutchk(0,6);
if ~isa(fun,'function_handle')
   error('FUN must be a function handle.');
end

if nargin>3
   validateattributes(options,{'nlopt.options'},{'scalar'},'nlopt.fminbnd','options');
   n = options.Dimension;
else
   n = numel(lb);
   try
      [~,~] = fun(x0);
      options = nlopt.options('gd_mlsl',n);
   catch
      options = nlopt.options('gn_direct',n);
   end
end

validateattributes(lb,{'double'},{'vector','nonnan'},'nlopt.fminbnd','x1');
validateattributes(ub,{'double'},{'vector','numel',n,'nonnan'},'nlopt.fminbnd','x2');
if any(lb(:)>=ub(:))
   error('Upper bounds must be strictly greater than the corresponding lower bounds');
end

% call the options.fmincon to perform the optimization
[varargout{1:nargout}] = options.fminbnd(fun,(lb+ub)/2,lb,ub);

end
