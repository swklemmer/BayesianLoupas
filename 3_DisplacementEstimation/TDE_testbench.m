addpath('../lib/DispEst/')

%%
% % Imaging parameters
% img_param = struct(...
%     'f_c',      7e6, ...    % [Hz] 
%     'bw',       0.5, ...    % [%]
%     't_s',      0, ...      % [s]
%     'z_max',    80, ...     % axial depth [lmbds]
%     'K',        0, ...      % number of kernels
%     'K_len',    3, ...      % kernel length [wvls]
%     'K_hop',    1, ...      % kernel hop [wvls]
%     'M',        0, ...      % kernel rf length [smpls]
%     'M_hop',    0, ...      % kernel rf hop [smpls]
%     'M_max',    0, ...      % axial depth [smpls]
%     'N',        3, ...      % ensemble length [frames]
%     'SNR',      10, ...     % [dB]
%     'SNR_rho',  0);         % [peak to peak]
% 
% img_param.t_s = 1 / (4 * img_param.f_c);
% img_param.SNR_rho = 10^(img_param.SNR/20);
% 
% % Method parameters
% met_param = struct(...
%     'lambda',   5e1, ...    % prior weight
%     'p',        1.05, ...   % norm deegree
%     'B',        5);         % vecinity size [kernels]

img_param = struct(...
    'fr_d', 1 ...         % frame rate decimation (1 --> 20 kHz)
    );

%%

% Load sonograms of different material properties
%ct_list = 0.25:0.25:3;  % shear wave speed [m/s]
ct_list = 1.75;
%N_exp = 20;              % nr. of experiment per sws
N_exp = 1;

for c_t = ct_list
    for n_i = 1:N_exp

        % Load results
        load(sprintf('../2_SonogramSimulation/results/d%d/ct%4.2f_%d.mat', ...
            50, c_t, n_i), ...
            'RcvData', 'IData', 'QData', 'img_x', 'img_z')

        % Beamform IQ data
        BFData = beamform_iq(IData{1}, QData{1}, img_param.fr_d);

        % Split line into kernels
%         [rf_data, img_param] = split_windows(img_param, rf_lines);


%         [rf_data, img_param] = split_windows(img_param, rf_lines);
%         
%         % Calculate starting solution using Loupas
%         iq_lines = demoulate_rf(img_param, rf_lines);
%         u_0 = loupas(img_param, iq_lines);
%         
%         % Maximize posterior probability
%         fun = @(u) -eval_posterior(img_param, met_param, rf_data, 'ack', u);
%         
%         opts = optimoptions('fmincon', ...
%             'Display', 'iter', ...
%             'MaxFunctionEvaluations', 5e4, ...
%             'OptimalityTolerance', 1e-4, ...
%             'StepTolerance', 1e-4);
%         u_hat = fmincon(fun, u_0, [], [], [], [],...
%             -0.5 * ones(img_param.K, 1),...
%             0.5 * ones(img_param.K, 1), [], opts);
    end
end