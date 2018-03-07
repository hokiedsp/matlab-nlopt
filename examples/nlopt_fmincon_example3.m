clear;

% Find the point where Rosenbrock's function is minimized within a circle, also subject to bound constraints.
fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% Look within the region $0 \le x(1) \le 0.5$, $0.2 \le x(2) \le 0.8$.

lb = [0,0.2];
ub = [0.5,0.8];

% There are no linear constraints, so set those arguments to [].
A = [];
b = [];
Aeq = [];
beq = [];

% Choose an initial point satisfying all the constraints.
x0 = [1/4,1/4];

% Solve the problem.
x = nlopt.fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@circlecon)



function c = circlecon(x)
% Look within the circle centered at [1/3,1/3] with radius 1/3. Copy the following code to a file on your MATLAB® path named circlecon.m.
c = (x(1)-1/3).^2 + (x(2)-1/3).^2 - (1/3).^2;
end
