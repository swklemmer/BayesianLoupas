% Experiment info
param_list = 0:5:60;

% Imaging parameters
img_param = struct(...
    'f_c',      7e6, ...    % [Hz] 
    'bw',       0.5, ...    % [%]
    't_s',      0, ...      % [s]
    'z_max',    10, ...     % [lmbds]
    'M',        0, ...      % [number of samples]
    'N',        3, ...      % [ensemble length]
    'SNR',      10, ...     % [dB]
    'SNR_rho',  1);         % [peak to peak]

img_param.t_s = 1 / (4 * img_param.f_c);
img_param.M = ceil(img_param.z_max / (img_param.f_c * img_param.t_s));

% Method parameters
met_param = struct(...
    'u_dim', -0.5:0.001:0.5, ... % [lmbds]
    'ack_a', 5e-2, ...           % [distribution width]
    'ncc_a', 5e-2);              % [distribution width]

% Simulation parameters
rng(420)
n_exp = 1e3; % [number of experiments]
u_true = rand([1, n_exp]) - 0.5; % [lmbds]

% Create gaussian RF Pulse 
[rf_pulse, t_cell] = create_pulse(img_param);

% Calculate SNR estimate
SNR = zeros(length(param_list), n_exp);

for p = 1:length(param_list)
    % Assign value to corresponing parameter
    [img_param, met_param] = update_parameter('SNR', param_list(p));

    for n = 1:n_exp % [lmbds]

        % Create RF lines
        rf_lines = create_rf_line(...
            img_param, rf_pulse, t_cell, u_true(n));

        % Estimate SNR
        SNR(p, n) = estimate_snr(rf_lines, img_param, met_param);
    end
end

%%
figure(1)
boxplot(20 * log10(SNR'), param_list)
hold on
plot(param_list, 5 * param_list, '--')
hold off
ylim([0, 60])
% figure(1)
% plot(param_list, 20 * log10(mean(SNR, 2)))
% hold on
% plot(param_list, param_list, '--')
% hold off
% ylim([param_list(1), param_list(end)])
% grid on
% 
% figure(2)
% plot(param_list, std(20 * log10(SNR), [], 2))
% ylim([0, param_list(end)])
% grid on