%clear all
addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')

% Experiment info
alg_list = {'ack'};

% Imaging parameters
img_param = struct(...
    'f_c',      7e6, ...    % [Hz] 
    'bw',       0.5, ...    % [%]
    't_s',      0, ...      % [s]
    'z_max',    8, ...      % [lmbds]
    'M',       -1, ...      % [number of samples]
    'N',        3, ...      % [ensemble length]
    'SNR',      5); ...    % [dB]

img_param.t_s = 1 / (4 * img_param.f_c);
img_param.M = ceil(img_param.z_max / (img_param.f_c * img_param.t_s));

% Method parameters
met_param = struct(...
    'u_dim', -.5:1e-3:.5, ... % [lmbds]
    'ack_a', 1e-3, ...           % [distribution width]
    'ncc_a', 2e-2, ...           % [distribution width]
    'qck_a', 1e-0);              % [distribution width]

% Simulation parameters
rng(19122023)
u_true = -0.3; % [lmbds]

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
    sum(p_xu) * diff(met_param.u_dim(1:2))
end
