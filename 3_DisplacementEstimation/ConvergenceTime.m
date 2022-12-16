addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')
clear all
graf = 0;

%% Simulation parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX', 'TW');
lambda = 1540 / (Trans.frequency * 1e6);

% Imaging parameters
img_param = struct(...
    'fr_n',     40, ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation (1 --> 20 kHz)
    'snr',      60, ...      % signal-to-noise-ratio [dB]
    'z_len',    8, ...      % axial kernel length [wvls]
    'z_hop',    2, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    2 ...       % temp. kernel length [wvls]
    );

% Method parameters
met_param = struct(...
    'alpha',    0, ...      % prior weight
    'p',        2, ...      % norm deegree
    'vcn_z',    3, ...      % axial vecinity size [kernels]
    'vcn_x',    1, ...      % late. vecinity size [kernels]
    'alg',      'ack' ...   % Likelihood function
    );

% Optimization parameters
opt_param = optimoptions('fmincon', ...
    'Display', 'none', ...
    'MaxFunctionEvaluations', 2e4, ...
    'MaxIterations', 3e1, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-6, ...
    'OutputFcn', @iteration_error ...
    );

%% Calculate convergence
snr_list = [5, 20, 60];

for snr = snr_list
    clear functions
    img_param.snr = snr;

    train_dir = sprintf('../resources/TrainData/SPW2/snr%d/', img_param.snr);
    train_files = char(dir(fullfile(train_dir, '*.mat')).name);
    
    % Pre-allocate convergence metrics
    err = zeros(size(train_files, 1), opt_param.MaxIterations, 3);
    err_j = zeros(opt_param.MaxIterations, 3);
    elapsed_t = zeros(size(train_files, 1), opt_param.MaxIterations);
    elapsed_t_j = zeros(opt_param.MaxIterations, 1);
    fun_eval = zeros(size(train_files, 1), opt_param.MaxIterations);
    fun_eval_j = zeros(opt_param.MaxIterations, 1);
    
    for j = 1:size(train_files, 1)
    
        fprintf('File %d\n', j);
    
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
    
        if j == 1
            % Pre-allocate estimations
            u_est = zeros([size(train_files, 1), size(u_0)]);
        end
    
        % Select optimal alpha
        met_param.alpha = optimal_alpha(img_param, met_param);
    
        % Select kernel corresponding to current frame
        fun = @(u) -eval_posterior_2D(met_param, RF_kern, u);
    
        % Maximize posterior probability for each frame
        u_hat = fmincon(fun, u_0, [], [], [], [],...
           -0.5 * ones(size(u_0)),...
            0.5 * ones(size(u_0)), [], opt_param);
    
        % Save final estimation
        u_est(j, :, :) = u_hat;
    
        % Save error and time metrics
        err(j, :, :) = err_j;
        elapsed_t(j, :) = elapsed_t_j;
        fun_eval(j, :) = fun_eval_j;
    end
    
    if ~graf
    save(sprintf('../resources/TimeData/%s/alpha%.0f_%ddB.mat', ...
         met_param.alg, met_param.p, img_param.snr), ...
        'u_est', 'err', 'elapsed_t', 'fun_eval', ...
        'img_param', 'met_param', 'opt_param')
    end
end
