fun = @(x)1+x(1,:)./(1+x(2,:)) - 3*x(1,:).*x(2,:) + x(2,:).*(1+x(1,:));

x = linspace(0,1,101);
y = linspace(0,2,101);
[X2,X1] = ndgrid(y,x);
F = zeros(size(X1));
F(:) = fun([X1(:) X2(:)].');
imagesc(x,y,F);

disp('nlopt.fmincon: Bound constraints only')
lb = [0,0];
ub = [1,2];
A = [];
b = [];
Aeq = [];
beq = [];
x0 = [0.5,1];
x = nlopt.fmincon(fun,x0,A,b,Aeq,beq,lb,ub)

disp('nlopt.fmincon: Try another initial x')
x0 = x0/5;
x = nlopt.fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
