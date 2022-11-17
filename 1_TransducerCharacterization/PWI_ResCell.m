addpath('../lib/Field_II/')
addpath('../lib/TransCharac/')
addpath('../lib/')
addpath('../../Vantage-4.7.6/Utilities')
rng(6942069)
field_init(0)
graf = 1;

%% Parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('L11-5v_PWI.mat', 'Trans', 'Parameters', 'TW');
Resource.Parameters = Parameters;

% Simulation parameters
img_param = struct(...
    'f_c',      TW(1).Parameters(1) * 1e6, ...  % Central frequency [Hz]
    'f_s',      250e6, ...                  % Sampling frequency [Hz]
    't_s',      1 / 250e6, ...              % Sampling frequency [Hz]
    'att',      0.3, ...                    % Attenuation [dB/cm/MHz]
    'n_ang',    3, ...                      % Steering angles
    'max_beta', 12, ...                     % Max angle [º]
    'kaiser',   3);                         % Kaiser window beta

%% Calculate resolution cell volume
param_list = [1:2:11];
val_list = zeros(length(param_list), 2);

for i = 1:length(param_list)
    img_param.n_ang = param_list(i);
    [v1, v2] = resolution_cell(img_param, Parameters, Trans, TW(1));
    val_list(i, :) = [v1, v2];
    fprintf("Res. Cell. Volume = %.4f mm^3\n", val_list(i, 1) * 1e9);
end

img_param.n_ang = 0;
[cv_foc, lr_foc] = resolution_cell(img_param, Parameters, Trans, TW(1));

%% Save results
save('results/steer_exp.mat', "param_list", 'val_list', 'img_param')

%% Show results
if graf
load('results/steer_exp.mat', "param_list", 'val_list', 'img_param')

fig = figure(1);
fig.Position = [0, 0, 500, 350];
plot(param_list, val_list(:, 1) * 1e9, '-*', 'LineWidth', 1)
ylim([0, 0.05])
yline(cv_foc * 1e9, '--')
xticks(param_list)
xlabel('# Angles', 'FontSize', 14)
ylabel('Cell volume [mm^3]', 'FontSize', 14)
title('-6dB Resolution Cell', 'FontSize', 16)
legend({'PWI', 'Focused'}, 'FontSize', 12)
grid on

fig = figure(2);
fig.Position = [500, 0, 500, 350];
plot(param_list, val_list(:, 2) * 1e3, '-*', 'LineWidth', 1)
ylim([0, 0.5])
yline(lr_foc * 1e3, '--')
xticks(param_list)
xlabel('# Angles', 'FontSize', 14)
ylabel('Lateral Beam Width [mm]', 'FontSize', 14)
title('Lateral resolution', 'FontSize', 16)
legend({'PWI', 'Focused'}, 'FontSize', 12)
grid on
end

%% Terminate program
field_end