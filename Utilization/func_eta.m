function [eta,deta] = func_eta(t)
    global a_ref f_ref
    eta_inf = 0.1;
    eta_0 = a_ref*8;
    % eta_inf = 10;
    % eta_0 = 15;
    a = 2.5;
    eta = (eta_0 - eta_inf) * exp(-a * t) + eta_inf;
    deta = -a * (eta_0 - eta_inf) * exp(-a * t);
end
