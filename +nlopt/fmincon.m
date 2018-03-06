function varargout = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,C,Ceq,options)
%nlopt.FMINCON Find constrained minimum of a multivariable function
%   nlopt.FMINCON attempts to solve problems of the form:
%
%    min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%     X                     C(X) <= 0, Ceq(X) = 0   (nonlinear constraints)
%                              LB <= X <= UB        (bounds)
%
%   nlopt.FMINCON defaults to LD_SLSQP algorithm if gradient is given by
%   the objective function or LN_COBYYA if no gradient is available. To use
%   a different algorithm, predefine a custom option:
%
%      OPTIONS = nlopt.options('YOUR_CHOICE',ndim)
%
%   See nlopt.options.getAlgorithms() for the complete list of algorithms,
%   and consult the official NLopt documentation for which algorithm
%   supports nonlinear constraints.
%
%   X = nlopt.FMINCON(FUN,X0,A,B) starts at X0 and finds a minimum X to the
%   function FUN, subject to the linear inequalities A*X <= B. FUN accepts
%   input X and returns a scalar function value F evaluated at X. X0 may be
%   a scalar, vector, or matrix.
%
%   X = nlopt.FMINCON(FUN,X0,A,B,Aeq,Beq) minimizes FUN subject to the linear
%   equalities Aeq*X = Beq as well as A*X <= B. (Set A=[] and B=[] if no
%   inequalities exist.)
%
%   X = nlopt.FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that a solution is found in
%   the range LB <= X <= UB. Use empty matrices for LB and UB
%   if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded below;
%   set UB(i) = Inf if X(i) is unbounded above.
%
%   X = nlopt.FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,C,Ceq) subjects the minimization
%   to the constraints defined in NONLCON. The function NONLCON accepts X
%   and returns the vectors C and Ceq, representing the nonlinear
%   inequalities and equalities respectively. nlopt.FMINCON minimizes FUN such
%   that C(X) <= 0 and Ceq(X) = 0. (Set LB = [] and/or UB = [] if no bounds
%   exist.)
%
%   X = nlopt.FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes with
%   the default optimization parameters replaced by values in OPTIONS, an
%   argument created with the OPTIMOPTIONS function. See OPTIMOPTIONS for
%   details. For a list of options accepted by nlopt.FMINCON refer to the
%   documentation.
%
%   [X,FVAL] = nlopt.FMINCON(FUN,X0,...) returns the value of the objective
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = nlopt.FMINCON(FUN,X0,...) returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are listed below. See the documentation
%   for a complete description.
%
%   All algorithms:
%     1  First order optimality conditions satisfied.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -2  No feasible point found.
%   Trust-region-reflective, interior-point, and sqp:
%     2  Change in X too small.
%   Trust-region-reflective:
%     3  Change in objective function too small.
%   Active-set only:
%     4  Computed search direction too small.
%     5  Predicted change in objective function too small.
%   Interior-point and sqp:
%    -3  Problem seems unbounded.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = nlopt.FMINCON(FUN,X0,...) returns a structure
%   OUTPUT with information such as total number of iterations, and final
%   objective function value. See the documentation for a complete list.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD] = nlopt.FMINCON(FUN,X0,...) returns the
%   value of the gradient of FUN at the solution X.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = nlopt.FMINCON(FUN,X0,...)
%   returns the value of the exact or approximate Hessian of the Lagrangian
%   at X.
%
%   Examples
%     FUN can be specified using @:
%        X = fmincon(@humps,...)
%     In this case, F = humps(X) returns the scalar function value F of
%     the HUMPS function evaluated at X.
%
%     FUN can also be an anonymous function:
%        X = fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0])
%     returns X = [0;0].
%
%   If FUN or NONLCON are parameterized, you can use anonymous functions to
%   capture the problem-dependent parameters. Suppose you want to minimize
%   the objective given in the function myfun, subject to the nonlinear
%   constraint mycon, where these two functions are parameterized by their
%   second argument a1 and a2, respectively. Here myfun and mycon are
%   MATLAB file functions such as
%
%        function f = myfun(x,a1)
%        f = x(1)^2 + a1*x(2)^2;
%
%        function [c,ceq] = mycon(x,a2)
%        c = a2/x(1) - x(2);
%        ceq = [];
%
%   To optimize for specific values of a1 and a2, first assign the values
%   to these two parameters. Then create two one-argument anonymous
%   functions that capture the values of a1 and a2, and call myfun and
%   mycon with two arguments. Finally, pass these anonymous functions to
%   nlopt.FMINCON:
%
%        a1 = 2; a2 = 1.5; % define parameters first
%        options = optimoptions('fmincon','Algorithm','interior-point'); % run interior-point algorithm
%        x = fmincon(@(x) myfun(x,a1),[1;2],[],[],[],[],[],[],@(x) mycon(x,a2),options)
%
%   See also nlopt.options, nlopt.options.getAlgorithms, nlopt.fminunc,
%            nlopt.fminbnd

narginchk(4,11);
nargoutchk(0,6);
if ~isa(fun,'function_handle')
   error('FUN must be a function handle.');
end

mode = false(1,4); % [Ab, Abeq lb ub C Ceq]

n = numel(x0);
if nargin > 2
   tf = [isempty(A) isempty(b)];
   if sum(tf)==1
      error('Linear inequality constraints must have both A and b.');
   end
   if tf(1)
      validateattributes(A,{'double'},{'2d','ncol',n,'finite'},'nlopt.fmincon','A');
      ncon = size(A,1);
      validateattributes(b,{'double'},{'column','nrows',ncon,'finite'},'nlopt.fmincon','b');
      mode(1) = true;
   end
   
   if nargin>4
      tf = [isempty(Aeq) nargin<6||isempty(beq)];
      if sum(tf)==1
         error('Linear equality constraints must have both Aeq and beq.');
      end
      if tf(1)
         validateattributes(Aeq,{'double'},{'2d','ncol',n,'finite'},'nlopt.fmincon','Aeq');
         ncon = size(Aeq,1);
         validateattributes(beq,{'double'},{'column','nrows',ncon,'finite'},'nlopt.fmincon','beq');
         mode(2) = true;
      end
      
      if nargin>6
         if ~isempty(lb)
            validateattributes(lb,{'double'},{'vector','numel',n,'nonnan'},'nlopt.fmincon','lb');
         end
         
         if nargin>7
            if ~isempty(ub)
               validateattributes(ub,{'double'},{'vector','numel',n,'nonnan'},'nlopt.fmincon','ub');
            end
            
            if nargin>8
               if ~isempty(C)
                  validateattributes(C,{'function_handle'},{},'nlopt.fmincon','C');
                  mode(3) = true;
               end
               if nargin>9
                  if ~isempty(Ceq)
                     validateattributes(Ceq,{'function_handle'},{},'nlopt.fmincon','Ceq');
                     mode(4) = true;
                  end
               end
            end
         else
            ub = [];
         end
      else
         lb = [];
      end
   end
end

if nargin>10
   validateattributes(options,{'nlopt.options'},{'scalar'},'nlopt.fmincon','options');
   validateattributes(x0,{'double'},{'vector','numel',options.Dimension,'real','finite'});
else
   validateattributes(x0,{'double'},{'vector','real','finite'});
   try
      [~,~] = fun(x0);
      options = nlopt.options('LD_SLSQP',n);
   catch
      options = nlopt.options('gn_isres',n);
   end
end

% combine them into lb, ub, con, mcon, coneq, mconeq
con = {};
mcon = {};
coneq = {};
mconeq = {};
if mode(1) % A&b given
   if isscalar(b)
      con{1} = @(x)A*x-b;
   else
      mcon{1} = @(x)A*x-b;
   end
end
if mode(2) %Aeq&beq
   if isscalar(beq)
      coneq{1} = @(x)Aeq*x-beq;
   else
      mconeq{1} = @(x)Aeq*x-beq;
   end
end
if mode(3) % C
   try
      val = C(x);
      if isscalar(val)
         con{end+1} = C;
      else
         mcon{end+1} = C;
      end
   catch
      error('Invalid nonlinear constraint function C(x)');
   end
end
if mode(4) % C
   try
      val = Ceq(x);
      if isscalar(val)
         coneq{end+1} = Ceq;
      else
         mconeq{end+1} = Ceq;
      end
   catch
      error('Invalid nonlinear equality constraint function Ceq(x)');
   end
end

% call the options.fmincon to perform the optimization
[varargout{1:nargout}] = options.fmincon(fun,x0,con,mcon,coneq,mconeq,lb,ub);
end
