clear;

anonrosen = @(x)(100*(x(2,:) - x(1,:).^2).^2 + (1-x(1,:)).^2);

x = linspace(-2,2,101);
y = linspace(-1,3,101);
[X2,X1] = ndgrid(y,x);
F = zeros(size(X1));
F(:) = anonrosen([X1(:) X2(:)].');
imagesc(x,y,F);

x0 = [-1;2];
% options = optimoptions(@fminunc,'Algorithm','quasi-newton');
[x, fval] = nlopt.fminunc(anonrosen,x0);

hold on
plot(x0(1),x0(2),'xk','MarkerFaceColor','k')
plot(x(1),x(2),'ok','MarkerFaceColor','k')
hold off
