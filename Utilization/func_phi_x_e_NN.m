function [phi,dphi]=func_phi_x_e_NN(x)
    x1 = x(1);
    x2 = x(2);
    e1 = x(3);
    e2 = x(4);
    % phi = [x1^2; x1*x2; x2^2; x1^2*x2; x1*x2^2; x1^2*x2^2];
    % dphi = [2*x1, 0; x2, x1; 0, 2*x2; 2*x1*x2, x1^2; x2^2, 2*x1*x2; 2*x1*x2^2, x1^2*x2];
    phi_x = [x1*e1; x1*e2; x2*e1; x2*e2; x1^2; x2^2]; % len = 6
    phi_e = [e1^2; e1*e2; e2^2]; % len = 3
    phi = [phi_x; phi_e];
    dphi_x = [x1 + e1 0; e2 x1; x2 e1; 0 x2 + e2; 2*x1 0; 0 2*x2];
    dphi_e = [2*e1, 0; e2, e1; 0, 2*e2];
    dphi = [dphi_x; dphi_e];