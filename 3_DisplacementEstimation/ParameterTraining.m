addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')
graf = 0;

%% Training parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX', 'TW');
lambda = 1540 / (Trans.frequency * 1e6);

% Imaging parameters
img_param = struct(...
    'fr_n',     40, ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation (1 --> 20 kHz)
    'snr',      5, ...     % signal-to-noise-ratio [dB]
    'z_len',    8, ...      % axial kernel length [wvls]
    'z_hop',    2, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    2 ...       % temp. kernel length [wvls]
    );

% Method parameters
met_param = struct(...
    'alpha',    3, ...     % prior weight
    'p',        1.05, ...      % norm deegree
    'vcn_z',    3, ...      % axial vecinity size [kernels]
    'vcn_x',    1 ...       % late. vecinity size [kernels]
    );

% Optimization parameters
opt_param = optimoptions('fmincon', ...
    'Display', 'none', ...
    'MaxFunctionEvaluations', 2e4, ...
    'MaxIterations', 3e1, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-6 ...
    );

%% Calculate errors

param_list = logspace(2.1, 4.5, 18);
%param_list = [0.25, 0.5, 1:20];
train_dir = sprintf('../resources/TrainData/SPW2/snr%d/', img_param.snr);
train_files = char(dir(fullfile(train_dir, '*.mat')).name);

% Pre-allocate error metrics
err = zeros(length(param_list), size(train_files, 1), 3);
err_0 = zeros(size(train_files, 1), 3);
a_gain = zeros(length(param_list), size(train_files, 1));

for j = 1:size(train_files, 1)

    % Load training data
    load([train_dir, train_files(j, : )], ...
        'RF_frames', 'I_frames', 'Q_frames', 'fr', 'c_t')
    
    % Split data into kernels
    [est_z, est_x, RF_kern, I_kern, Q_kern] = ...
        split_kernels(img_param, PData, RF_frames, I_frames, Q_frames);

    % Load true displacement
    [u_tru, ~]= interp_real_u(est_x, est_z, fr, PData, lambda, c_t);

    % Calculate starting solution using Loupas and it's error metrics
    u_0 = loupas_3D(I_kern, Q_kern);
    err_0(j, :) = error_metrics(u_0, u_tru, 1);

    if j == 1
        % Pre-allocate estimations
        u_est = zeros([length(param_list), size(train_files, 1), size(u_0)]);
    end

    for i = 1:length(param_list); tic();

        % Update parameter
        met_param.alpha = param_list(i);

        % Select kernel corresponding to current frame
        fun = @(u) -eval_posterior_2D(met_param, RF_kern, 'ack', u);

        % Maximize posterior probability for each frame
        u_hat = fmincon(fun, u_0, [], [], [], [],...
           -0.5 * ones(size(u_0)),...
            0.5 * ones(size(u_0)), [], opt_param);

        % Find best gain
        a_gain(i, j) = norm(u_hat(:) - u_tru(:), 2) / norm(u_hat(:),  2);
    
        % Calculate error metrics (without gain)
        err(i, j, :) = error_metrics(u_hat, u_tru, 1);

        fprintf('i = %d | c_t =%5.2f | t =%4.0f\n', i, c_t, toc());

        % Save estimation
        u_est(i, j, :, :) = u_hat;

        if graf
            % Show estimations versus ground truth
            compare_frame(est_x, est_z, u_0, u_hat, u_tru)
            pause()
        end
    end
end

if ~graf
save(sprintf('../resources/ErrorMetrics/alpha4_%ddB.mat', img_param.snr), ...
    'u_est', 'err', 'err_0', 'a_gain', ...
    'param_list', 'img_param', 'met_param', 'opt_param')
end
