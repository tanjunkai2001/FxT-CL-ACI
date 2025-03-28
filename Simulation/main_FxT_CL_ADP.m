close all
clear
clc

cd('../')
addpath('Utilization')
addpath('Data')
addpath('Simulation')
%% 系统parameters
global len_x len_u len_d len_W1 len_W2 len_t a_ref f_ref
len_x = 2;
len_u = 2;
len_d = 1;
len = 6;
% len = 3;
len_W1 = len;
len_W2 = len;
len_t = 20;
T = 0.0005;
% x = [0.5; 0.5];
a_ref = 0.5; f_ref = 2;
[xref, dxref] = func_reference(0);
% x = xref + a_ref * 5 * rand(len_x, 1);
x = xref + a_ref * 5 * ones(len_x, 1);
e = x - xref;

global model len_theta dW_theta_ a_theta Gamma_theta
model = 0; % moel-free, SysID
% model = 1; % model-based
len_theta = 4;
W_theta = 0.5 * ones(len_theta, 1);
W_theta_ = [-1 1 -0.5 -0.5]';
a_theta = 6;
% Gamma_theta = diag([0.5, 0.5, 1.1, 1.1]);
Gamma_theta = diag([0.5, 0.5, 0.5, 0.5, ]);
dW_theta_ = zeros(len_theta, 10);


%% 控制器parameters
global gamma Q R a_c1 a_c2 a_u dW_c1_ dW_c2_ dW_a1_ dW_a2_ N temp umax U V H2Hinf Bellman
H2Hinf = 1;
umax = 5;
gamma = 5;
Q = 1 * eye(len_x);
% R = 5 * eye(len_u);
R = 20 * eye(len_u);

a_c1 = 0.5;
a_c2 = 0.5;
a_u = 1;

W0 = 2;
W_c1 = W0 * ones(len_W1, 1);
W_c2 = W0 * ones(len_W2, 1);
W_a1 = W0 * ones(len_W1, 1);
W_a2 = W0 * ones(len_W2, 1);
% W_c1 = W0 * rand(len_W1, 1);
% W_c2 = W0 * rand(len_W2, 1);
% W_a1 = W0 * rand(len_W1, 1);
% W_a2 = W0 * rand(len_W2, 1);

y = 0.5 * ones(len_x, 1);
N = 1;
dW_c1_ = zeros(len_W1, N);
dW_c2_ = zeros(len_W2, N);
dW_a1_ = zeros(len_W1, N);
dW_a2_ = zeros(len_W2, N);
temp = 1;

U = [];
V = [];
Bellman = [];

%% ode
% X0 = [x; e; W_c1; W_c2; W_a1; W_a2; y];
X0 = [x; e; W_c1; W_c2; W_a1; W_a2; W_theta; y];
tspan = 0:T:len_t;
% options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
% [t, X] = ode23(@func_ode, tspan, X0, options);
[t, X] = ode23(@func_ode, tspan, X0);
% [t, X] = RK1(@func_ode, [0 len_t], X0, T);
% [t, X] = RK4(@func_ode, [0 len_t], X0, T);
% X = X';
% t = t';

x = X(:, 1:len_x);
e = X(:, len_x + 1:2 * len_x);
W_c1 = X(:, 2 * len_x + 1:2 * len_x + len_W1);
W_c2 = X(:, 2 * len_x + len_W1 + 1:2 * len_x + len_W1 + len_W2);
W_a1 = X(:, 2 * len_x + len_W1 + len_W2 + 1:2 * len_x + len_W1 + 2 * len_W1);
W_a2 = X(:, 2 * len_x + len_W1 + len_W2 + len_W1 + 1:2 * len_x + len_W1 + 2 * len_W1 + len_W2);
% y = X(:, 2 * len_x + len_W1 + len_W2 + 2 * len_W1 + 1:2 * len_x + len_W1 + len_W2 + 2 * len_W1 + len_x);
W_theta = X(:, 2 * len_x + len_W1 + len_W2 + 2 * len_W1 + 1:2 * len_x + len_W1 + len_W2 + 2 * len_W1 + len_theta);
y = X(:, 2 * len_x + len_W1 + len_W2 + 2 * len_W1 + len_theta + 1:2 * len_x + len_W1 + len_W2 + 2 * len_W1 + len_theta + len_x);

xref = x - e;

save('Data/data_case1.mat');

%% plot fig
% plot_fig_FxT_CL_ADP

%% ode function
function dX = func_ode(t, X_)
    % System
    global gamma Q R a_c1 a_c2 a_u len_x len_W1 len_W2 len_t len_u len_d U V
    % PPC
    global alpha_l alpha_h
    % Learning
    global dW_c1_ dW_c2_ dW_a1_ N temp umax H2Hinf Bellman
    % SysID
    global model len_theta dW_theta_ a_theta Gamma_theta
    x = X_(1:len_x);
    e = X_(len_x + 1:2 * len_x);
    W_c1 = X_(2 * len_x + 1:2 * len_x + len_W1);
    W_c2 = X_(2 * len_x + len_W1 + 1:2 * len_x + len_W1 + len_W2);
    W_a1 = X_(2 * len_x + len_W1 + len_W2 + 1:2 * len_x + len_W1 + 2 * len_W1);
    W_a2 = X_(2 * len_x + len_W1 + len_W2 + len_W1 + 1:2 * len_x + len_W1 + 2 * len_W1 + len_W2);
    W_theta = X_(2 * len_x + len_W1 + len_W2 + 2 * len_W1 + 1:2 * len_x + len_W1 + len_W2 + 2 * len_W1 + len_theta);
    y = X_(2 * len_x + len_W1 + len_W2 + 2 * len_W1 + len_theta + 1:2 * len_x + len_W1 + len_W2 + 2 * len_W1 + len_theta + len_x);

    %% Controller
    t
    % [xref, dxref] = func_reference(t);
    % e = x - xref;
    xref = x - e;
    dxref = [-1 1; -2 1] * xref;

    x1 = x(1);
    x2 = x(2);

    if model == 1
        % f = [x2 - 2 * x1; -x2 - 0.5 * x1 + 0.25 * x2 * cos(x1) ^ 2 + 0.25 * x2 * sin(2 * x1 + 1) ^ 2];
        % g = [0; cos(x1)];
        % k = [0; sin(2 * x1)];
        f = [-x1 + x2; -0.5 * x1 - 0.5 * x2 * (1 - (cos(2 * x1 + 1)) ^ 2)];
        f_ = f;
    else
        W_theta_ = [-1 1 -0.5 -0.5]';
        f = func_fx(x) * W_theta_;
        f_ = func_fx(x) * W_theta;
    end

    g = diag([sin(2 * x1 + 1) + 2, cos(2 * x1) + 2]);
    k = 1 * ones(2, 1);
    k = 1 * ones(2, 1);

    %% NN basis: simple NN, StaF NN (exp), StaF NN (modified)
    % [phi1, dphi1] = Staf_Basis(epsilon);
    % [phi2, dphi2] = Staf_Basis(epsilon);
    % [phi1, dphi1] = Staf_Basis_modified(epsilon);
    % [phi2, dphi2] = Staf_Basis_modified(epsilon);
    % [phi1, dphi1] = func_phi(e);
    % [phi2, dphi2] = func_phi(e);
    [phi1, dphi1] = func_phi_6NN(e);
    [phi2, dphi2] = func_phi_6NN(e);
    % [phi1, dphi1] = func_phi_8NN(e);
    % [phi2, dphi2] = func_phi_8NN(e);
    % [phi1, dphi1] = func_phi_x_e_NN([x; e]);
    % [phi2, dphi2] = func_phi_x_e_NN([x; e]);

    % u = W_a1' * phi;
    % u = -0.5 * R \ (Psi*g)' * dphi1' * W_a1;
    % v = 0.5 / gamma ^ 2 * k' * dphi2' * W_a2;
    u = -umax * tanh(0.5 / umax * R \ (g)' * dphi1' * W_a1);
    % u = umax * tanh(W_a1' * phi1 / umax);
    v = umax * tanh(0.5 / umax * k' * dphi2' * W_a2);

    % if t < 3
    %     signal = 0.0 * exp(-t / 2) * tanh(t / 5) * [(sin(t) + cos(2 * t) + sin(3 * t) + cos(4 * t) + sin(5 * t) +cos(6 * t) + sin(7 * t) +cos(8 * t) + sin(9 * t)); (cos(t) + sin(2 * t) + cos(3 * t) + sin(4 * t) + cos(5 * t) + sin(6 * t) + cos(7 * t) + sin(8 * t) + cos(9 * t))];
    % else
    %     signal = [0; 0];
    % end

    % v = umax * tanh(signal(1) / umax);
    % u = u + 2 * pinv(g) * (dxref - f_);
    dx = f + g * u + k * v;
    de = dx - dxref;
    de_ = f_ + g * u + k * v - dxref;

    U = [U; [u', t]];
    V = [V; [v', t]];

    dy = zeros(len_x, 1);

    %% stackelberg
    df = [1 * W_theta(1), 1 * W_theta(2); 1 * W_theta(3) + W_theta(4) * 2 * x2 * (-sin(2 * x1)), W_theta(4) * (cos(2 * x1) + 2)];
    dg1 = [2 * cos(2 * x1 + 1) 0; 0 0];
    dg2 = [0 0; 0 0];
    dk = [0 0; 0 0];
    dk1 = [0; 0];
    dk2 = [0; 0];
    % dk = [0 0; 2 * cos(2 * x1), 0];
    % dk1 = [0; 2 * cos(2 * x1)];
    % dk2 = [0; 0];

    % du = dphi1' * W_a1;
    % dphi1 = [2*x1, 0; x2, x1; 0, 2*x2; 2*x1*x2, x1^2; x2^2, 2*x1*x2; 2*x1*x2^2, x1^2*x2];
    % ddphi1 = [2, 0; 0, 0; 0, 2; 0, 0; 0, 0; 0, 0];

    du1 = [0; 0];
    du2 = [0; 0];
    dv_dV2 = -0.5 / gamma ^ 2 * k';
    dV1 = dphi1' * W_c1;
    dV2 = dphi1' * W_c2;
    % size(du)
    dlambda2 =- ((df + dg1 * u(1) + dg2 * u(2) + dk * v + g(:, 1) * du1' + g(:, 2) * du2')' * dV2 - 2 * Q * e - 2 * du1 * R(1) * u(1) - 2 * du2 * R(2) * u(2));
    dy =- ((k * dv_dV2)' * dV1 - (df + dg1 * u(1) + dg2 * u(2) + dk * v + g(:, 1) * du1' + g(:, 2) * du2') * y - (y(1) * dv_dV2' * dk1' * dV2 + y(2) * dv_dV2' * dk2' * dV2));

    %% update
    gamma1 = dphi1 * de_;
    gamma2 = dphi2 * de_;

    if H2Hinf == 0
        r1 = e' * Q * e + u' * R * u - gamma ^ 2 * v' * v;
        r2 = r1;
    else
        r1 = e' * Q * e + u' * R * u;
        r2 = r1 - gamma ^ 2 * v' * v;
    end

    % dy = [x' * Q * x + u' * R * u; r2];
    % r2 = gamma ^ 2 * v' * v-r1;
    % r2 = gamma ^ 2 * v' * v - r1;
    delta1 = W_c1' * gamma1 + r1 + 1 * y' * dlambda2;
    delta2 = W_c2' * gamma2 + r2;
    delta3 = (f - f_);
    delta4 = (W_a1 - W_c1);
    delta5 = (W_a2 - W_c2);


    gamma1_FxT = 0.8;
    gamma2_FxT = 1.2;

    delta1_FxT = (abs(delta1)) ^ gamma1_FxT * sign(delta1) + (abs(delta1)) ^ gamma2_FxT * sign(delta1);
    delta2_FxT = (abs(delta2)) ^ gamma1_FxT * sign(delta2) + (abs(delta2)) ^ gamma2_FxT * sign(delta2);
    delta3_FxT = (abs(delta3)) .^ gamma1_FxT .* sign(delta3) + (abs(delta3)) .^ gamma2_FxT .* sign(delta3);
    delta4_FxT = (abs(delta4)) .^ gamma1_FxT .* sign(delta4) + (abs(delta4)) .^ gamma2_FxT .* sign(delta4);
    delta5_FxT = (abs(delta5)) .^ gamma1_FxT .* sign(delta5) + (abs(delta5)) .^ gamma2_FxT .* sign(delta5);

    
    Bellman = [Bellman; [r1, r2, t]];
    dW_c1_current = -a_c1 * gamma1 * delta1_FxT / (gamma1' * gamma1 + 1) ^ 2;
    dW_c2_current = -a_c2 * gamma2 * delta2_FxT / (gamma2' * gamma2 + 1) ^ 2;
    % dW_a1_current = -a_u * (W_a1 - W_c1);
    % dW_a2_current = -a_u * (W_a2 - W_c2);
    dW_a1_current = -a_u * delta4_FxT;
    dW_a2_current = -a_u * delta5_FxT;
    dW_theta_current = Gamma_theta * a_theta * func_fx(x)' * delta3_FxT;

    dW_c1_(:, temp) = dW_c1_current;
    dW_c2_(:, temp) = dW_c2_current;
    dW_a1_(:, temp) = dW_a1_current;
    dW_a2_(:, temp) = dW_a2_current;
    dW_theta_(:, temp) = dW_theta_current;
    temp = mod(temp, N) + 1;
    dW_c1 = sum(dW_c1_, 2) + 9 * dW_c1_current;
    dW_c2 = sum(dW_c2_, 2) + 9 * dW_c2_current;
    dW_a1 = sum(dW_a1_, 2) + 9 * dW_a1_current;
    dW_a2 = sum(dW_a2_, 2) + 9 * dW_a2_current;
    dW_theta = sum(dW_theta_, 2) + 0 * dW_theta_current;

    dX = [dx; de; dW_c1; dW_c2; dW_a1; dW_a2; dW_theta; dy];
end
