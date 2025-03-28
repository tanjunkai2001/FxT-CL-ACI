function f_basis = func_fx_2(x)
    x1 = x(1);
    x2 = x(2);
    f_basis = [x1; x2; x2 * (cos(2 * x1) + 2)];
    % f_basis = [x1 x2 0 0
    %            0 0 x1 x2 * (1 - (cos(2 * x1 + 1)) ^ 2)];
    % f_basis = [x1 x2 0 0
    %            0 0 x1 x2 * (cos(2 * x1) + 2)];
