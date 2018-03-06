clear; close all; drawnow;

opt = nlopt.options('Algorithm', 'LD_MMA', 1,'StepRelativeTolerance',1e-4,'ConstraintTolerance',
opt.lower_bounds = [-inf, 0]
opt.min_objective = @myfunc
opt.fc = { (@(x) myconstraint(x,2,0)), (@(x) myconstraint(x,-1,1)) }
opt.fc_tol = [1e-8, 1e-8];
opt.xtol_rel = 1e-4

[xopt, fmin, retcode] = nlopt_optimize(opt, [1.234 5.678])


function [val, gradient] = myfunc(x)
    val = sqrt(x(2));
    if (nargout > 1)
        gradient = [0, 0.5 / val];
    end
    end

    function [val, gradient] = myconstraint(x,a,b)
    val = (a*x(1) + b)^3 - x(2);
    if (nargout > 1)
        gradient = [3*a*(a*x(1) + b)^2, -1];
    end
    end
