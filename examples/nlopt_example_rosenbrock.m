clear; close all; drawnow

anonrosen = @(x)(100*(x(2,:) - x(1,:).^2).^2 + (1-x(1,:)).^2);
x0 = [-1;2];

x = linspace(-2,2,101);
y = linspace(-1,3,101);
[X2,X1] = ndgrid(y,x);
F = zeros(size(X1));
F(:) = anonrosen([X1(:) X2(:)].');
imagesc(x,y,F);
line(x0(1),x0(2),'Marker','x','LineStyle','none','MarkerFaceColor','w','Color','w');
h = line(x0(1),x0(2),'Marker','.','Color','k');

options = nlopt.options('ln_neldermead',2,'OutputFun',@(x,optimvalues,state)outfun(x,optimvalues,state,h));
[x, fval, flag, output] = nlopt.fminunc(anonrosen,x0,options);

h1 = line(x0(1),x0(2),'Marker','.','Color','r','LineStyle','--');
options = nlopt.options('ld_lbfgs',2,'OutputFun',@(x,optimvalues,state)outfun(x,optimvalues,state,h1));
[x, fval, flag, output] = nlopt.fminunc(@rosenbrockwithgrad,x0,options);

line(x(1),x(2),'Marker','o','LineStyle','none','Color','w','MarkerFaceColor','w')

legend([h h1],'Nelder-Mead', 'Limited-memory BFGS')


function stop = outfun(x,~,~,h)
  set(h,'XData',[h.XData x(1)],'YData',[h.YData x(2)]);
  stop = false;
end

function [f,g] = rosenbrockwithgrad(x)
% Calculate objective f
f = 100*(x(2,:) - x(1,:).^2).^2 + (1-x(1,:)).^2;

if nargout > 1 % gradient required
    g = [-400*(x(2,:)-x(1,:).^2).*x(1,:)-2*(1-x(1,:));
        200*(x(2,:)-x(1,:).^2)];
end
end
