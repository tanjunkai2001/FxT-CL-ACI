function [phi,dphi]=func_phi(x)
    x1 = x(1);
    x2 = x(2);
    phi = [x1^2; x1*x2; x2^2];
    dphi = [2*x1, 0; x2, x1; 0, 2*x2];
