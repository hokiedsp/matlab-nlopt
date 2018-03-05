clear;

anonrosen = @(x)(100*(x(2,:) - x(1,:).^2).^2 + (1-x(1,:)).^2);
x0 = [-1;2];

x = linspace(-2,2,101);
y = linspace(-1,3,101);
[X2,X1] = ndgrid(y,x);
F = zeros(size(X1));
F(:) = anonrosen([X1(:) X2(:)].');
imagesc(x,y,F);
line(x0(1),x0(2),'Marker','x','LineStyle','none','MarkerFaceColor','k','Color','k');
h = line(x0(1),x0(2),'Marker','.','Color','k');

options = nlopt.options('ln_neldermead',2,'OutputFun',@outfun);
[x, fval, flag, output] = nlopt.fminunc(anonrosen,x0,options);

line(x(1),x(2),'Marker','o','LineStyle','none','Color','k','MarkerFaceColor','k')

function stop = outfun(x,optimvalues,state)
  set(h,'XData',[h.XData x(1)],'YData',[h.YData x(2)]);
  stop = false;
end
