addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')
clear all
graf = 0;

%% Simulation parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_PWI.mat', ...
    'PData', 'Parameters', 'Trans', 'TX', 'TW', 'Receive');
lambda = 1540 / (TW(1).Parameters(1) * 1e6);

% Imaging parameters
img_p = struct(...
    'f_c',      TW(1).Parameters(1) * 1e6, ...  % central frequency [Hz]
    'f_s',      Receive(1).ADCRate * 1e6, ...   % sampling frequency [Hz]
    't_s',      -1, ...                         % sampling frequency [Hz]
    'snr',      -1 ...                         % signal-to-noise-ratio [dB]
    );

img_p.t_s = 1 / img_p.f_s;

% Estimation parameters
est_p = struct(...
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    1, ...      % axial kernel hop    [wvls]
    'x_len',    .5, ...      % late. kernel length [wvls]
    'x_hop',    .5, ...      % late. kernel hop    [wvls]
    't_len',    2, ...      % temp. kernel length [frames]
    'crs_mref', 0, ...      % coarse est. moving reference frame
    'fin_mref', 1, ...      % fine est. moving reference frame
    'crs_win',  1, ...      % coarse est. windowing
    'fin_win',  0 ...      % fine est. windowing
   );

% Method parameters
met_p = struct(...
    'alpha',   -1, ...      % prior weight
    'p',        1, ...      % norm deegree
    'vcn_z',    3, ...      % axial vecinity size [kernels]
    'vcn_x',    1, ...      % late. vecinity size [kernels]
    'alg',      'ack' ...   % Likelihood function
    );

% Optimization parameters
opt_p = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...
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
    img_p.snr = snr;

    train_dir = sprintf('../resources/TrainData/SPW/snr%d/', img_p.snr);
    train_files = char(dir(fullfile(train_dir, '*.mat')).name);
    
    % Pre-allocate convergence metrics
    err = zeros(size(train_files, 1), opt_p.MaxIterations, 3);
    err_j = zeros(opt_p.MaxIterations, 3);
    elapsed_t = zeros(size(train_files, 1), opt_p.MaxIterations);
    elapsed_t_j = zeros(opt_p.MaxIterations, 1);
    fun_eval = zeros(size(train_files, 1), opt_p.MaxIterations);
    fun_eval_j = zeros(opt_p.MaxIterations, 1);
    
    for j = 1:size(train_files, 1)
    
        fprintf('File %d\n', j);
    
        % Load training data
        load(strtrim([train_dir, train_files(j, : )]), ...
            'RF_frames', 'I_frames', 'Q_frames', 'fr', 'c_t');
        
        % Split data into kernels
        [est_z, est_x, RF_kern, I_kern, Q_kern] = ...
            split_kernels(est_p, PData, RF_frames, I_frames, Q_frames);
    
        % Load true displacement
        [u_tru, ~]= interp_real_u(est_x, est_z, fr, PData, lambda, c_t);
    
        % Calculate starting solution using Loupas and it's error metrics
        u_0 = loupas_3D(est_p, I_kern, Q_kern);
    
        if j == 1
            % Pre-allocate estimations
            u_est = zeros([size(train_files, 1), size(u_0)]);
        end
    
        % Select optimal alpha
        met_p.alpha = optimal_alpha(img_p, met_p);
    
        % Select kernel corresponding to current frame
        fun = @(u) -eval_posterior_2D(img_p, met_p, RF_kern, u);
    
        % Maximize posterior probability for each frame
        u_hat = fmincon(fun, u_0, [], [], [], [],...
           -0.5 * ones(size(u_0)),...
            0.5 * ones(size(u_0)), [], opt_p);
    
        % Save final estimation
        u_est(j, :, :) = u_hat;
    
        % Save error and time metrics
        err(j, :, :) = err_j;
        elapsed_t(j, :) = elapsed_t_j;
        fun_eval(j, :) = fun_eval_j;
    end
    
    if ~graf
    save(sprintf('../resources/TimeData/%s/alpha%.0f_%ddB.mat', ...
         met_p.alg, met_p.p, img_p.snr), ...
        'u_est', 'err', 'elapsed_t', 'fun_eval', ...
        'img_p', 'met_p', 'opt_p')
    end
end
