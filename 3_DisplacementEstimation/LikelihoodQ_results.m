clear all


%% Kernel axial length

% Experiment info
exp_param = 'z_max';
img_snr = 20;

% Load results
load(sprintf('results/%s_%ddB_4.mat', exp_param, img_snr),...
    'res_q', 'param_list', 'alg_list');

fig = figure();
fig.Position = [0, 100, 350, 230];
plot(param_list, res_q)
ylim([0, 1])
title('Effect of kernel length on likelihood quality', 'FontSize', 12)
xlabel('Axial kernel length [\lambda]', 'FontSize', 12)
ylabel('Median Probability', 'FontSize', 12)
legend(alg_list, 'FontSize', 12, 'Position', [0.78 0.21 0.05 0.05])
grid on

%% SNR

% Experiment info
exp_param = 'alpha8';
img_snr = 20;

% Load results
load(sprintf('results/%s_%ddB.mat', exp_param, img_snr),...
    'res_q', 'param_list', 'alg_list');

fig = figure();
fig.Position = [0, 100, 350, 230];
semilogx(param_list, res_q)
title('Effect of \alpha on likelihood quality', 'FontSize', 12)
ylim([0, 1])
xlabel('Parameter \alpha', 'FontSize', 12)
ylabel('Median Probability', 'FontSize', 12)
legend(alg_list, 'FontSize', 12)
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
