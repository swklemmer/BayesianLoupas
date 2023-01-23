clear all


%% Kernel axial length

% Experiment info
exp_param = 'z_max';
img_snr = 20;

% Load results
load(sprintf('results/%s_%ddB_3.mat', exp_param, img_snr),...
    'res_q', 'param_list', 'alg_list');

fig = figure();
fig.Position = [0, 100, 260, 230];
plot(param_list, res_q(:, 1))
hold on
plot(param_list, res_q(:, 2), '--')
hold off
ylim([0, 1])
title('Effect of kernel length on likelihood quality', 'FontSize', 12)
%title(' ', 'FontSize', 12)
xlabel('Axial kernel length [\lambda]', 'FontSize', 12)
ylabel('Median Probability', 'FontSize', 12)
legend({'ACK', 'NCC'}, 'FontSize', 12, 'Location', 'southeast')
grid on

%% SNR

% Experiment info
exp_param = 'alpha8';
img_snr = 20;

% Load results
load(sprintf('results/%s_%ddB.mat', exp_param, img_snr),...
    'res_q', 'param_list', 'alg_list');

fig = figure();
fig.Position = [0, 100, 250, 230];
semilogx(param_list, res_q(:, 1))
hold on
semilogx(param_list, res_q(:, 2), '--')
hold off
title('Effect of \beta on likelihood quality', 'FontSize', 12)
%title(' ', 'FontSize', 12)
ylim([0, 1])
xlabel('Parameter \beta', 'FontSize', 12)
ylabel('Median Probability', 'FontSize', 12)
legend({'ACK', 'NCC'}, 'FontSize', 12)
grid on

%% Computation Time

fig = figure();
fig.Position = [500, 100, 300, 200];
semilogx(param_list, res_t)
ylim([0, 7])
title(sprintf('Algorithm Time (%s)', exp_param), 'Interpreter','none')
xlabel(exp_param, 'Interpreter','none')
ylabel('Computation Time [ms]')
legend(alg_list)
grid on
