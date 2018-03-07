clear;

fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

disp('nlopt.fmincon: Linear inequality constraint only')
x0 = [-1,2];
A = [1,2];
b = 1;
x = nlopt.fmincon(fun,x0,A,b)

disp('nlopt.fmincon: Linear inequality and equality constraints')
x0 = [0.5,0];
A = [1,2];
b = 1;
Aeq = [2,1];
beq = 1;
x = nlopt.fmincon(fun,x0,A,b,Aeq,beq)
