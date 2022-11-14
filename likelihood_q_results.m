clear all

% Experiment info
exp_param = 'alpha_5';

% Load results
load(sprintf('QualityResults/%s.mat', exp_param),...
    'res_q', 'res_t', 'param_list', 'alg_list', 'img_param');

fig = figure();
fig.Position = [0, 100, 300, 200];
semilogx(param_list, res_q)
ylim([0, 1])
title(sprintf('Likelihood Quality (%s)', exp_param))
xlabel(exp_param)
ylabel('Median Probability')
legend(alg_list)
grid on

fig = figure();
fig.Position = [500, 100, 300, 200];
semilogx(param_list, res_t)
ylim([0, 7])
title(sprintf('Algorithm Time (%s)', exp_param))
xlabel(exp_param)
ylabel('Computation Time [ms]')
legend(alg_list)
grid on
