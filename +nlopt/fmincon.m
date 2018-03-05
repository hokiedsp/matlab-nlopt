function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
%FMINCON finds a constrained minimum of a function of several variables.
%   FMINCON attempts to solve problems of the form:
%    min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%     X                     C(X) <= 0, Ceq(X) = 0   (nonlinear constraints)
%                              LB <= X <= UB        (bounds)
%    
%   FMINCON implements four different algorithms: interior point, SQP,
%   active set, and trust region reflective. Choose one via the option
%   Algorithm: for instance, to choose SQP, set OPTIONS =
%   optimoptions('fmincon','Algorithm','sqp'), and then pass OPTIONS to
%   FMINCON.
%                                                           
%   X = FMINCON(FUN,X0,A,B) starts at X0 and finds a minimum X to the 
%   function FUN, subject to the linear inequalities A*X <= B. FUN accepts 
%   input X and returns a scalar function value F evaluated at X. X0 may be
%   a scalar, vector, or matrix. 
%
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq) minimizes FUN subject to the linear 
%   equalities Aeq*X = Beq as well as A*X <= B. (Set A=[] and B=[] if no 
%   inequalities exist.)
%
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that a solution is found in 
%   the range LB <= X <= UB. Use empty matrices for LB and UB
%   if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; 
%   set UB(i) = Inf if X(i) is unbounded above.
%
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) subjects the minimization
%   to the constraints defined in NONLCON. The function NONLCON accepts X 
%   and returns the vectors C and Ceq, representing the nonlinear 
%   inequalities and equalities respectively. FMINCON minimizes FUN such 
%   that C(X) <= 0 and Ceq(X) = 0. (Set LB = [] and/or UB = [] if no bounds
%   exist.)
%
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes with
%   the default optimization parameters replaced by values in OPTIONS, an
%   argument created with the OPTIMOPTIONS function. See OPTIMOPTIONS for
%   details. For a list of options accepted by FMINCON refer to the
%   documentation.
%  
%   X = FMINCON(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the linear inequality constraints in PROBLEM.Aineq
%   and PROBLEM.bineq, the linear equality constraints in PROBLEM.Aeq and
%   PROBLEM.beq, the lower bounds in PROBLEM.lb, the upper bounds in 
%   PROBLEM.ub, the nonlinear constraint function in PROBLEM.nonlcon, the
%   options structure in PROBLEM.options, and solver name 'fmincon' in
%   PROBLEM.solver. Use this syntax to solve at the command line a problem 
%   exported from OPTIMTOOL. 
%
%   [X,FVAL] = FMINCON(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = FMINCON(FUN,X0,...) returns an EXITFLAG that
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
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINCON(FUN,X0,...) returns a structure 
%   OUTPUT with information such as total number of iterations, and final 
%   objective function value. See the documentation for a complete list.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = FMINCON(FUN,X0,...) returns the 
%   Lagrange multipliers at the solution X: LAMBDA.lower for LB, 
%   LAMBDA.upper for UB, LAMBDA.ineqlin is for the linear inequalities, 
%   LAMBDA.eqlin is for the linear equalities, LAMBDA.ineqnonlin is for the
%   nonlinear inequalities, and LAMBDA.eqnonlin is for the nonlinear 
%   equalities.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD] = FMINCON(FUN,X0,...) returns the 
%   value of the gradient of FUN at the solution X.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = FMINCON(FUN,X0,...) 
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
%   FMINCON:
%
%        a1 = 2; a2 = 1.5; % define parameters first
%        options = optimoptions('fmincon','Algorithm','interior-point'); % run interior-point algorithm
%        x = fmincon(@(x) myfun(x,a1),[1;2],[],[],[],[],[],[],@(x) mycon(x,a2),options)
%
%   See also OPTIMOPTIONS, OPTIMTOOL, FMINUNC, FMINBND, FMINSEARCH, @, FUNCTION_HANDLE.

%   Copyright 1990-2015 The MathWorks, Inc.
