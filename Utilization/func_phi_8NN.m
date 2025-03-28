function [phi,dphi]=func_phi_8NN(x)
    x1 = x(1);
    x2 = x(2);
    phi = [x1^2; x1*x2; x2^2; x1^2*x2; x1*x2^2; x1^2*x2^2; x1^3; x2^3];
    dphi = [2*x1, 0; x2, x1; 0, 2*x2; 2*x1*x2, x1^2; x2^2, 2*x1*x2; 2*x1*x2^2, x1^2*x2; 3*x1^2, 0; 0, 3*x2^2];
