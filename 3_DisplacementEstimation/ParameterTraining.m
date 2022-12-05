addpath('../lib/DispEst/')
addpath('../lib/SonoSim/BMS_aux/')
graf = 0;

%% Training parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX', 'TW');
lambda = 1540 / (TW(1).Parameters(1) * 1e6);

% Imaging parameters
img_param = struct(...
    'fr_n',     40, ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation (1 --> 20 kHz)
    'snr',      20, ...     % signal-to-noise-ratio [dB]
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    5, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    3 ...       % temp. kernel length [wvls]
    );

% Method parameters
met_param = struct(...
    'lambda',   1e2, ...    % prior weight
    'p',        2, ...      % norm deegree
    'vcn_z',    3, ...      % axial vecinity size [kernels]
    'vcn_x',    3 ...       % late. vecinity size [kernels]
    );         

%% Optimize parameters
param_list = [1e0, 5e0, 1e1, 5e1, 1e2, 5e2];
train_dir = sprintf('../resources/TrainData/snr%d/', img_param.snr);
train_files = char(dir(fullfile(train_dir, '*.mat')).name);

% Pre-allocate error metrics
e_bias = zeros(length(param_list), size(train_files, 1));
e_var = zeros(length(param_list), size(train_files, 1));
e_rmse = zeros(length(param_list), size(train_files, 1));

for i = 1:length(param_list)
    met_param.lambda = param_list(i);
    for j = 1:size(train_files, 1)
    
        % Load training data
        load([train_dir, train_files(j, : )], ...
            'BM_frames', 'I_frames', 'Q_frames', 'u_rea', 'rea_x', 'rea_z')
    
        % Split data into kernels
        [est_z, est_x, RF_kern, I_kern, Q_kern] = ...
            split_kernels(img_param, PData, BM_frames, I_frames, Q_frames);
    
        % Calculate starting solution using Loupas
        u_0 = loupas_3D(I_kern, Q_kern);
    
        % Maximize posterior probability for each frame
        opts = optimoptions('fmincon', ...
            'Display', 'iter-detailed', ...
            'MaxFunctionEvaluations', 2e4, ...
            'MaxIterations', 3e1, ...
            'OptimalityTolerance', 1e-4, ...
            'StepTolerance', 1e-4, ...
            'TypicalX', u_0);
    
        % Select kernel corresponding to current frame
        fun = ...
            @(u) -eval_posterior_2D(met_param, squeeze(RF_kern), 'ack', u);
    
        % Solve bayesian regression
        u_hat = fmincon(fun, u_0, [], [], [], [],...
           -0.5 * ones(size(u_0)),...
            0.5 * ones(size(u_0)), [], opts);
    
        % Interpolate in space dimention
        [x_grid, z_grid] = meshgrid(est_x * lambda, est_z * lambda);
        u_tru = interp2(rea_z, rea_x, u_rea, z_grid(:), x_grid(:), ...
            'linear', 0);
        u_tru = reshape(u_tru, size(x_grid)) / lambda;
    
        % Calculate error metrics
        e_bias(i, j) = mean(u_hat - u_tru, 'all');
        e_var(i, j) = var(u_hat - u_tru, [], 'all');
        e_rmse(i, j) = sqrt(e_var(i, j) + e_bias(i, j)^2);
    
        if graf
            % Show estimations versus ground truth
            compare_frame(est_x, est_z, u_0, u_hat, u_tru)
            pause()
        end
    end
end
