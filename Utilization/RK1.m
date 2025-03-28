%% Forward euler method for solving differential equations
function [t, x] = Simul_Dyn(ufunc, tspan, x0, h)
    % input arguments: 'function name'--ufun，'time starting and ending points'--tspan,'initial state'--x0，'step size'--h
    % output arguments: 'time'--t，'state'--x，'internal state'--eta, 'lyapunov function'--V, 'augmented function'--W

    if nargin < 4 % the number of input arguments
        h = 0.0001;
    end

    if size(tspan) == [1, 2]
        t0 = tspan(1);
        tn = tspan(2);
    else
        error(message('MATLAB:Euler:WrongDimensionOfTspan'));
    end

    n = floor((tn - t0) / h);
    t(1) = t0;
    x(:, 1) = x0;

    for i = 1:n
        t(i + 1) = t(i) + h;
        [dx] = ufunc(t(i), x(:, i));
        x(:, i + 1) = x(:, i) + h * dx;
    end
