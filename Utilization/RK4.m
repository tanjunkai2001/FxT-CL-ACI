%龙格库塔四阶解微分方程
function [t, x] = RK4(f, tspan, x0, h)
    %   f:微分方程组表达式
    %   x0:方程组变量初值
    %   h:设定步长
    %   t0-t1:求解区间

    if nargin < 4 % the number of input arguments
        h = 0.001;
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

        k1 = f(t(i), x(:, i));
        k2 = f(t(i) + h / 2, x(:, i) + h / 2 * k1);
        k3 = f(t(i) + h / 2, x(:, i) + h / 2 * k2);
        k4 = f(t(i) + h, x(:, i) + h * k3);

        dx = (1/6) * (k1 + 2 * k2 + 2 * k3 + k4);
        
        x(:, i + 1) = x(:, i) + h * dx;
    end

end
