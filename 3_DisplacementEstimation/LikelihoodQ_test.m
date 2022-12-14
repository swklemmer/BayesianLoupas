clear all
addpath('../lib/SonoSim')
addpath('../lib/DispEst')

% Experiment info
alg_list = {'ack', 'ncc'};
exp_param = 'z_max';
param_list = 2:15;

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX', 'TW', 'Receive');
lambda = 1540 / (Trans.frequency * 1e6);

% Imaging parameters
img_param = struct(...
    'f_c',      Trans.frequency * 1e6, ...      % [Hz] 
    'bw',       diff(Trans.Bandwidth) / Trans.frequency, ...    % [%]
    't_s',      1e-6 / Receive(1).decimSampleRate, ...          % [s]
    'z_max',    8, ...                          % [lmbds]
    'M',        0, ...                          % [number of samples]
    'N',        3, ...                          % [ensemble length]
    'SNR',      20 ...                          % [dB]
    );

img_param.M = ceil(img_param.z_max / (img_param.f_c * img_param.t_s));

% Method parameters
met_param = struct(...
    'u_dim', -0.5:0.001:0.5, ... % [lmbds]
    'ack_a', 4.6e-3, ...           % [distribution width]
    'ncc_a', 1.8e-3);              % [distribution width]

% Simulation parameters
rng(420)
n_exp = 1e3; % [number of experiments]
u_true = max(min(randn([1, n_exp]) / 3 , 0.5), -0.5); % [lmbds]
eval_q = zeros(length(param_list), length(alg_list), n_exp);
eval_t = zeros(length(param_list), length(alg_list), n_exp);
res_q = zeros(length(param_list), length(alg_list));
res_t = zeros(length(param_list), length(alg_list));

% Create gaussian RF Pulse 
[rf_pulse, t_cell] = create_pulse(img_param);

for p = 1:length(param_list)

    % Assign value to corresponing parameter
    [img_param, met_param] = update_parameter(exp_param, param_list(p));

    for alg = 1:length(alg_list)

        % Raise parameter change flag
        param_flag = 1;

        for n = 1:n_exp % [lmbds]
    
            % Create RF lines
            rf_lines = create_rf_line(...
                img_param, rf_pulse, t_cell, u_true(n));

            % Calculate likelihood function using selected algorithm
            [p_xu, elapsed_t] = likelihood(...
                img_param, met_param, alg_list{alg}, rf_lines);

            % Integrate probability over +- lambda/100
            eval_q(p, alg, n) = likelihood_quality(...
                p_xu, met_param.u_dim, u_true(n), 1/100);
        
            % Save processing time
            eval_t(p, alg, n) = elapsed_t;

        end

        % Process results
        res_q(p, alg) = median(eval_q(p, alg, :));
        res_t(p, alg) = mean(eval_t(p, alg, :)) * 1e3;
        fprintf('%s @ p=%d: Q =%7.3f, (t =%5.1f ms)\n', ...
            alg_list{alg}, p, res_q(p, alg), res_t(p, alg))
    end
end

% Save results
save(sprintf('results/%s_%ddB_3.mat', exp_param, img_param.SNR),...
    'res_q', 'res_t', 'param_list', 'alg_list', 'img_param');
