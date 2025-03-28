close all
addpath('utilization');
addpath('Data');
addpath('Figures');
addpath('../ZoomPlot-MATLAB');

load('../Data/data_case1.mat')

markerindices = 1:len_t / T / 10:length(t) - 1;
markersize = 7;

savefig_flag = 0;

zoom_flag = 1;
tightfig_flag = 0;
linewidth = 1.5;
fontname = 'Times New Roman';
fontsize = 12;
figuresize = [100, 100, 500, 450];
figuresize2 = [100, 100, 550, 450]

path_fig = '../Figures/';
path_pwd = pwd;
path_now = pwd;

max(e(:, 1))
max(e(:, 2))
min(e(:, 1))
min(e(:, 2))
y(end,:)


%% weights 1
figure
set(gcf,'Color','w');
colorSet = lines(6);
set(gcf, 'Position', figuresize);

t1 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% Plot W_c1
nexttile
hold on
linestyles = {'-', '--', ':', '-.', '-', '--'};  % Distinct line styles
markers = {'o', 's', 'v', '^', 'd', 'p'}; % Different markers for each line

for i = 1:6
    plot(t, W_c1(:,i), 'Color', colorSet(i,:), 'LineWidth', 1.5,...
         'LineStyle', linestyles{i},...
         'Marker', markers{i}, 'MarkerIndices', markerindices,...
         'MarkerSize', 6)
end
ylabel('$\hat{W}_{c1}$','Interpreter','latex','FontSize',12)
legend('$\hat{W}_{c1}(1)$','$\hat{W}_{c1}(2)$','$\hat{W}_{c1}(3)$','$\hat{W}_{c1}(4)$','$\hat{W}_{c1}(5)$','$\hat{W}_{c1}(6)$',...
       'Interpreter','latex','Location','eastoutside',...
       'NumColumns',1,'FontSize',10,'box','off')
ylim([min(W_c1(end,:))-0.5, max(W_c1(end,:))+0.5])
grid on; box on
set(gca,'GridLineStyle','--','FontSize',10)

% Plot W_a1
nexttile
hold on
for i = 1:6
    plot(t, W_a1(:,i), 'Color', colorSet(i,:), 'LineWidth', 1.5,...
         'LineStyle', linestyles{i},...
         'Marker', markers{i}, 'MarkerIndices', markerindices,...
         'MarkerSize', 6)
end
ylabel('$\hat{W}_{a1}$','Interpreter','latex','FontSize',12)
legend('$\hat{W}_{a1}(1)$','$\hat{W}_{a1}(2)$','$\hat{W}_{a1}(3)$','$\hat{W}_{a1}(4)$','$\hat{W}_{a1}(5)$','$\hat{W}_{a1}(6)$',...
       'Interpreter','latex','Location','eastoutside',...
       'NumColumns',1,'FontSize',10,'box','off')
ylim([min(W_a1(end,:))-0.5, max(W_a1(end,:))+0.5])
grid on; box on
set(gca,'GridLineStyle','--','FontSize',10)

% % Plot W_c2
% nexttile
% hold on
% for i = 1:6
%     plot(t, W_c2(:,i), 'Color', colorSet(i,:), 'LineWidth', 1.5,...
%          'LineStyle', linestyles{i},...
%          'Marker', markers{i}, 'MarkerIndices', markerindices,...
%          'MarkerSize', 6)
% end
% ylabel('$\hat{W}_{c2}$','Interpreter','latex','FontSize',12)
% legend('$\hat{W}_{c2}(1)$','$\hat{W}_{c2}(2)$','$\hat{W}_{c2}(3)$','$\hat{W}_{c2}(4)$','$\hat{W}_{c2}(5)$','$\hat{W}_{c2}(6)$',...
%        'Interpreter','latex','Location','northoutside',...
%        'NumColumns',6,'FontSize',6,'box','off')
% ylim([min(W_c2(end,:))-0.5, max(W_c2(end,:))+0.5])
% grid on; box on
% set(gca,'GridLineStyle','--','FontSize',10)

% Plot W_theta with improved colors
nexttile
hold on
colorSet_theta = [0 0.4470 0.7410;  % blue
                 0.8500 0.3250 0.0980;  % orange
                 0.9290 0.6940 0.1250;  % yellow
                 0.4940 0.1840 0.5560;  % purple
                 0.4660 0.6740 0.1880;  % green
                 0.3010 0.7450 0.9330]; % light blue
for i = 1:len_theta
    plot(t, W_theta(:,i), 'LineWidth', 1.5, 'Color', colorSet_theta(i,:),...
         'LineStyle', linestyles{i},...
         'Marker', markers{i}, 'MarkerIndices', markerindices, 'MarkerSize', 6)
    yline(W_theta_(i), ':', 'LineWidth', 1.2, 'Color', colorSet_theta(i,:), 'Alpha', 0.3)
end
xlabel('Time [s]','FontSize',12,'FontName',fontname)
ylabel('$\hat{W}_{\theta}$','Interpreter','latex','FontSize',12)
legend('$\hat{W}_{\theta}(1)$', '', '$\hat{W}_{\theta}(2)$', '' ,'$\hat{W}_{\theta}(3)$', '', '$\hat{W}_{\theta}(4)$',...
       'Interpreter','latex','Location','eastoutside',...
       'NumColumns',1,'FontSize',10,'box','off')
grid on; box on
set(gca,'GridLineStyle','--','FontSize',10)

% 获取当前图窗大小
fig = gcf;
fig.Units = 'inches';
figPosition = fig.Position;
width = figPosition(3);
height = figPosition(4);

% 设置保存选项
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);

% 保存为 PDF 文件
% path_fig = '../Figure/';
% path_now = pwd;
cd(path_fig)
if savefig_flag == 1
    print(gcf, 'weights1.pdf', '-dpdf');
end
cd(path_pwd);


%% weights 2
figure
set(gcf,'Color','w');
colorSet = lines(6);
set(gcf, 'Position', figuresize);

t1 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% Plot W_c2
nexttile
hold on
linestyles = {'-', '--', ':', '-.', '-', '--'};  % Distinct line styles
markers = {'o', 's', 'v', '^', 'd', 'p'}; % Different markers for each line

for i = 1:6
    plot(t, W_c2(:,i), 'Color', colorSet(i,:), 'LineWidth', 1.5,...
         'LineStyle', linestyles{i},...
         'Marker', markers{i}, 'MarkerIndices', markerindices,...
         'MarkerSize', 6)
end
ylabel('$\hat{W}_{c2}$','Interpreter','latex','FontSize',12)
legend('$\hat{W}_{c2}(1)$','$\hat{W}_{c2}(2)$','$\hat{W}_{c2}(3)$','$\hat{W}_{c2}(4)$','$\hat{W}_{c2}(5)$','$\hat{W}_{c2}(6)$',...
       'Interpreter','latex','Location','eastoutside',...
       'NumColumns',1,'FontSize',10,'box','off')
ylim([min(W_c2(end,:))-0.5, max(W_c2(end,:))+0.5])
grid on; box on
set(gca,'GridLineStyle','--','FontSize',10)

% Plot W_a2
nexttile
hold on
for i = 1:6
    plot(t, W_a2(:,i), 'Color', colorSet(i,:), 'LineWidth', 1.5,...
         'LineStyle', linestyles{i},...
         'Marker', markers{i}, 'MarkerIndices', markerindices,...
         'MarkerSize', 6)
end
ylabel('$\hat{W}_{a2}$','Interpreter','latex','FontSize',12)
legend('$\hat{W}_{a2}(1)$','$\hat{W}_{a2}(2)$','$\hat{W}_{a2}(3)$','$\hat{W}_{a2}(4)$','$\hat{W}_{a2}(5)$','$\hat{W}_{a2}(6)$',...
       'Interpreter','latex','Location','eastoutside',...
       'NumColumns',1,'FontSize',10,'box','off')
ylim([min(W_a2(end,:))-0.5, max(W_a2(end,:))+0.5])
grid on; box on
set(gca,'GridLineStyle','--','FontSize',10)

% % Plot W_c2
% nexttile
% hold on
% for i = 1:6
%     plot(t, W_c2(:,i), 'Color', colorSet(i,:), 'LineWidth', 1.5,...
%          'LineStyle', linestyles{i},...
%          'Marker', markers{i}, 'MarkerIndices', markerindices,...
%          'MarkerSize', 6)
% end
% ylabel('$\hat{W}_{c2}$','Interpreter','latex','FontSize',12)
% legend('$\hat{W}_{c2}(1)$','$\hat{W}_{c2}(2)$','$\hat{W}_{c2}(3)$','$\hat{W}_{c2}(4)$','$\hat{W}_{c2}(5)$','$\hat{W}_{c2}(6)$',...
%        'Interpreter','latex','Location','northoutside',...
%        'NumColumns',6,'FontSize',6,'box','off')
% ylim([min(W_c2(end,:))-0.5, max(W_c2(end,:))+0.5])
% grid on; box on
% set(gca,'GridLineStyle','--','FontSize',10)

% Plot W_theta with improved colors
nexttile
hold on
colorSet_theta = [0 0.4470 0.7410;  % blue
                 0.8500 0.3250 0.0980;  % orange
                 0.9290 0.6940 0.1250;  % yellow
                 0.4940 0.1840 0.5560;  % purple
                 0.4660 0.6740 0.1880;  % green
                 0.3010 0.7450 0.9330]; % light blue
for i = 1:len_theta
    plot(t, W_theta(:,i), 'LineWidth', 1.5, 'Color', colorSet_theta(i,:),...
         'LineStyle', linestyles{i},...
         'Marker', markers{i}, 'MarkerIndices', markerindices, 'MarkerSize', 6)
    yline(W_theta_(i), ':', 'LineWidth', 1.2, 'Color', colorSet_theta(i,:), 'Alpha', 0.3)
end
xlabel('Time [s]','FontSize',12,'FontName',fontname)
ylabel('$\hat{W}_{\theta}$','Interpreter','latex','FontSize',12)
legend('$\hat{W}_{\theta}(1)$', '', '$\hat{W}_{\theta}(2)$', '' ,'$\hat{W}_{\theta}(3)$', '', '$\hat{W}_{\theta}(4)$',...
       'Interpreter','latex','Location','eastoutside',...
       'NumColumns',1,'FontSize',10,'box','off')
grid on; box on
set(gca,'GridLineStyle','--','FontSize',10)

% 获取当前图窗大小
fig = gcf;
fig.Units = 'inches';
figPosition = fig.Position;
width = figPosition(3);
height = figPosition(4);

% 设置保存选项
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);

% 保存为 PDF 文件
% path_fig = '../Figure/';
% path_now = pwd;
cd(path_fig)
if savefig_flag == 1
    print(gcf, 'weights2.pdf', '-dpdf');
end
cd(path_pwd);




%% Control input
figure
set(gcf,'Color','w');
% Use a more professional color scheme
colorSet = [0 0.4470 0.7410;     % Deep blue
           0.8500 0.3250 0.0980]; % Dark orange
set(gcf, 'Position', figuresize);

t_ = U(:,3);
u_ = reshape(U(:,1:len_x), 2, length(t_))';
v_ = reshape(V(:,1:len_d), 1, length(t_))';

t1 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% Control input 1
nexttile
hold on
% Plot control bounds first (in background)
fill([t_; flipud(t_)], [ones(size(t_))*umax; ones(size(t_))*(-umax)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
% Plot control signal 1
plot(t_, u_(:, 1), 'Color', colorSet(1,:), 'LineWidth', 1.5, ...
    'LineStyle', '-')
% Plot bounds
plot(t_, ones(size(t_))*umax, 'r:', 'LineWidth', 1.2)
plot(t_, ones(size(t_))*(-umax), 'r:', 'LineWidth', 1.2, 'HandleVisibility', 'off')

ylim([-umax*1.1, umax*1.1])
ylabel('Control input $u_1$','Interpreter','latex', 'FontSize', 12)
legend('', '$u(1)$', 'Bound', 'Interpreter','latex', 'Location','northeast', 'FontSize', 12, 'box', 'off', 'NumColumns', 2)
% legend('', '$u(1)$', 'Bound', 'Interpreter','latex', 'Location','eastoutside', 'FontSize', 12, 'box', 'off', 'NumColumns', 1)
grid on
set(gca,'GridLineStyle','--','FontSize',10)
box on

% Control input 2
nexttile
hold on
% Plot control bounds
fill([t_; flipud(t_)], [ones(size(t_))*umax; ones(size(t_))*(-umax)], ...
    [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3)
% Plot control signal 2
plot(t_, u_(:, 2), 'Color', colorSet(2,:), 'LineWidth', 1.5, ...
    'LineStyle', '-')
% Plot bounds
plot(t_, ones(size(t_))*umax, 'r:', 'LineWidth', 1.2)
plot(t_, ones(size(t_))*(-umax), 'r:', 'LineWidth', 1.2, 'HandleVisibility', 'off')

ylim([-umax*1.1, umax*1.1])
ylabel('Control input $u_2$','Interpreter','latex', 'FontSize', 12)
legend('', '$u(2)$', 'Bound', 'Interpreter','latex', 'Location','northeast', 'FontSize', 12, 'box', 'off', 'NumColumns', 2)
% legend('', '$u(2)$', 'Bound', 'Interpreter','latex', 'Location','eastoutside', 'FontSize', 12, 'box', 'off', 'NumColumns', 1)
grid on
set(gca,'GridLineStyle','--','FontSize',10)
box on

% Disturbance input subplot
nexttile
hold on
plot(t_, v_(:, 1), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, ...
    'LineStyle', '-')
xlabel('Time [s]','FontSize',12,'FontName',fontname)
ylabel('Disturbance $v$','Interpreter','latex', 'FontSize', 12)
legend('$v$', 'Interpreter','latex', 'Location','best', 'FontSize', 12, 'box', 'off')
% legend('$v$', 'Interpreter','latex', 'Location','eastoutside', 'FontSize', 12, 'box', 'off')
grid on
set(gca,'GridLineStyle','--','FontSize',10)
box on


% 获取当前图窗大小
fig = gcf;
fig.Units = 'inches';
figPosition = fig.Position;
width = figPosition(3);
height = figPosition(4);

% 设置保存选项
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);

% 保存为 PDF 文件
% path_fig = '../Figure/';
% path_now = pwd;
cd(path_fig)
if savefig_flag == 1
    print(gcf, 'control.pdf', '-dpdf');
end
cd(path_pwd);




%% y and Bellman error
figure
set(gcf,'Color','w');
set(gcf, 'Position', figuresize2);

% Professional color scheme
colorSet = [0 0.4470 0.7410;      % Deep blue
           0.8500 0.3250 0.0980;   % Dark orange
           0.4660 0.6740 0.1880;   % Green
           0.4940 0.1840 0.5560];  % Purple

% Output subplot
subplot(2, 1, 1)
plot(t, y(:,1), 'Color', colorSet(1,:), 'LineWidth', 1.5, 'LineStyle', '-')
hold on
plot(t, y(:,2), 'Color', colorSet(2,:), 'LineWidth', 1.5, 'LineStyle', '--')
grid on
set(gca, 'GridLineStyle', '--', 'FontSize', 10)
ylabel('Output $y$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$y_1$', '$y_2$', 'Location', 'best', 'NumColumns', 2,...
       'FontSize', 15, 'Interpreter', 'latex', 'box', 'off')
box on

if savefig_flag && zoom_flag
    zp = BaseZoom();
    zp.run;
end

% Bellman error subplot
subplot(2, 1, 2)
t_ = Bellman(:,3);
plot(t_, Bellman(:,1), 'Color', colorSet(3,:), 'LineWidth', 1.5, 'LineStyle', '-')
hold on
plot(t_, Bellman(:,2), 'Color', colorSet(4,:), 'LineWidth', 1.5, 'LineStyle', '--')
grid on
set(gca, 'GridLineStyle', '--', 'FontSize', 10)
xlabel('Time [s]', 'FontSize', 12, 'FontName', fontname)
ylabel('Bellman error', 'FontSize', 12, 'FontName', fontname)
legend('$H_2$ Bellman error', '$H_{\infty}$ Bellman error',...
       'Interpreter', 'latex', 'Location', 'best',...
       'NumColumns', 2, 'FontSize', 12, 'box', 'off')
box on

if savefig_flag && zoom_flag
    zp = BaseZoom();
    zp.run;
end

% 获取当前图窗大小
fig = gcf;
fig.Units = 'inches';
figPosition = fig.Position;
width = figPosition(3);
height = figPosition(4);

% 设置保存选项
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);

% 保存为 PDF 文件
% path_fig = '../Figure/';
% path_now = pwd;
cd(path_fig)
if savefig_flag == 1
    print(gcf, 'y_bellman.pdf', '-dpdf');
end
cd(path_pwd);



%%
figure
set(gcf,'Color','w');
set(gcf, 'Position', figuresize2);

% Professional color scheme
colorSet = [0 0.4470 0.7410;      % Deep blue
           0.8500 0.3250 0.0980];  % Dark orange

% Subplot for x1
subplot(2, 1, 1)
plot(t, xref(:,1), 'Color', colorSet(2,:), 'LineWidth', 1.5, 'LineStyle', '-')
hold on
plot(t, x(:,1), 'Color', colorSet(1,:), 'LineWidth', 1.5, 'LineStyle', '--')
grid on
set(gca, 'GridLineStyle', '--', 'FontSize', 10)
ylabel('State $x_1$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$x_{1,d}$', '$x_1$', 'Location', 'best', 'NumColumns', 2,...
       'FontSize', 15, 'Interpreter', 'latex', 'box', 'off')
box on

% if savefig_flag && zoom_flag
%     zp = BaseZoom();
%     zp.run;
% end

% Subplot for x2
subplot(2, 1, 2)
plot(t, xref(:,2), 'Color', colorSet(2,:), 'LineWidth', 1.5, 'LineStyle', '-')
hold on
plot(t, x(:,2), 'Color', colorSet(1,:), 'LineWidth', 1.5, 'LineStyle', '--')
grid on
set(gca, 'GridLineStyle', '--', 'FontSize', 10)
xlabel('Time [s]', 'FontSize', 12, 'FontName', fontname)
ylabel('State $x_2$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$x_{2,d}$', '$x_2$', 'Location', 'best', 'NumColumns', 2,...
       'FontSize', 15, 'Interpreter', 'latex', 'box', 'off')
box on

% if savefig_flag && zoom_flag
%     zp = BaseZoom();
%     zp.run;
% end

% 获取当前图窗大小
fig = gcf;
fig.Units = 'inches';
figPosition = fig.Position;
width = figPosition(3);
height = figPosition(4);

% 设置保存选项
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);

% 保存为 PDF 文件
% path_fig = '../Figure/';
% path_now = pwd;
cd(path_fig)
if savefig_flag == 1
    print(gcf, 'x1_x2.pdf', '-dpdf');
end
cd(path_pwd);




%% tracking error
figure
set(gcf,'Color','w');
set(gcf, 'Position', figuresize2);

% Professional color scheme
colorSet = [0 0.4470 0.7410;      % Deep blue
           0.8500 0.3250 0.0980];  % Dark orange

% Subplot for tracking error e1
subplot(2, 1, 1)
plot(t, e(:,1), 'Color', colorSet(1,:), 'LineWidth', 1.5, 'LineStyle', '-')
grid on
set(gca, 'GridLineStyle', '--', 'FontSize', 10)
ylabel('Tracking error $e_1$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$e_1$', 'Location', 'best', 'FontSize', 15, 'Interpreter', 'latex', 'box', 'off')
box on

if savefig_flag && zoom_flag
    zp = BaseZoom();
    zp.run;
end

% Subplot for tracking error e2
subplot(2, 1, 2)
plot(t, e(:,2), 'Color', colorSet(2,:), 'LineWidth', 1.5, 'LineStyle', '-')
grid on
set(gca, 'GridLineStyle', '--', 'FontSize', 10)
xlabel('Time [s]', 'FontSize', 12, 'FontName', fontname)
ylabel('Tracking error $e_2$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$e_2$', 'Location', 'best', 'FontSize', 15, 'Interpreter', 'latex', 'box', 'off')
box on

if savefig_flag && zoom_flag
    zp = BaseZoom();
    zp.run;
end

% Save figure
fig = gcf;
fig.Units = 'inches';
figPosition = fig.Position;
width = figPosition(3);
height = figPosition(4);

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);

cd(path_fig)
if savefig_flag == 1
    print(gcf, 'e1_e2.pdf', '-dpdf');
end
cd(path_pwd);