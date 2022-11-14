clear all
addpath('lib/')
addpath('lib/fraccircshift/')

% Experiment info
alg_list = {'ack'};

% Imaging parameters
img_param = struct(...
    'f_c',      7e6, ...    % [Hz] 
    'bw',       0.5, ...    % [%]
    't_s',      0, ...      % [s]
    'z_max',    8, ...     % [lmbds]
    'M',        0, ...      % [number of samples]
    'N',        3, ...      % [ensemble length]
    'SNR',      10, ...     % [dB]
    'SNR_rho',  1);         % [peak to peak]

img_param.t_s = 1 / (4 * img_param.f_c);
img_param.M = ceil(img_param.z_max / (img_param.f_c * img_param.t_s));
img_param.SNR_rho = 10^(img_param.SNR/20);

% Method parameters
met_param = struct(...
    'u_dim', -0.5:1e-3:0.5, ... % [lmbds]
    'ack_a', 5e-3, ...           % [distribution width]
    'ncc_a', 2e-2);              % [distribution width]

% Simulation parameters
rng(420)
u_true = -0.2; % [lmbds]

% Create gaussian RF Pulse 
[rf_pulse, t_cell] = create_pulse(img_param);

% Create RF lines
rf_lines = create_rf_line(img_param, rf_pulse, t_cell, u_true);

for alg = 1:length(alg_list)

    % Raise parameter change flag
    param_flag = 1;

    % Calculate likelihood function using selected algorithm
    [p_xu, elapsed_t] = likelihood(...
        img_param, met_param, alg_list{alg}, rf_lines, u_true);

    % Integrate 
    sum(p_xu) * (met_param.u_dim(2) - met_param.u_dim(1));
end
