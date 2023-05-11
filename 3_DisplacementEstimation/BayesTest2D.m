% This scripts estimates the displacements of a single US frame collection
% using a bayesian approach. Data is fetched from the "TrainFrames" folder.

addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')

%% Simulation parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_PWI.mat', ...
    'PData', 'Parameters', 'Trans', 'TX', 'TW', 'Receive');
lambda = 1540 / (TW(1).Parameters(1) * 1e6);

% Imaging parameters
img_p = struct(...
    'f_c',      TW(1).Parameters(1) * 1e6, ...  % central frequency [Hz]
    'f_s',      TW(1).Parameters(1) / PData.PDelta(3) * 1e6, ...   % sampling frequency [Hz]
    't_s',      -1, ...                         % sampling frequency [Hz]
    'snr',      60 ...                         % signal-to-noise-ratio [dB]
    );

img_p.t_s = 1 / img_p.f_s;

% Estimation parameters
est_p = struct(...
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    2, ...     % axial kernel hop    [wvls]
    'x_len',    .5, ...    % late. kernel length [wvls]
    'x_hop',    2, ...     % late. kernel hop    [wvls]
    't_len',    2, ...     % temp. kernel length [frames]
    'fin_mref', 1, ...     % fine est. moving reference frame
    'fin_win',  1 ...      % fine est. windowing
   );

% Method parameters
met_p = struct(...
    'alg',      'sqck', ...   % Likelihood function
    'p',        2, ...      % norm deegree
    'alpha',    10^(-1), ...      % prior weight
    'vcn_z',    3, ...      % axial vecinity size [kernels]
    'vcn_x',    1 ...      % late. vecinity size [kernels]
    );

% Optimization parameters
opt_p = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...
    'MaxFunctionEvaluations', 2e4, ...
    'MaxIterations', 10, ...
    'OptimalityTolerance', 1e-2, ...
    'StepTolerance', 1e-5 ...
    );

%% Run estimation

% Load frames
c_t = 1.25; % [m/s]
it = 1;
load(sprintf('../resources/TrainData/SPW/snr%d/ct%.2f_%d.mat', ...
    img_p.snr, c_t, it), ...
    'RF_frames', 'I_frames', 'Q_frames', 'fr', 'c_t')

% Split data into kernels
[est_z, est_x, RF_kern, I_kern, Q_kern] = ...
    split_kernels(est_p, PData, RF_frames, I_frames, Q_frames);

% Load ground truth
load(sprintf('../resources/EstData/ground_truth/%.2f.mat', c_t), 'u_mean');

% Calculate starting solution and its error metrics
u_0 = loupas_3D(est_p, I_kern, Q_kern);
err_0 = error_metrics(u_0, u_mean, 1);

% Pre-compute best fast frequency columns
max_fc = pre_comp_SQCK(img_p, u_0, RF_kern);
param_flag = 1;

% Declare loss function
fun = @(u) -eval_posterior_2D(img_p, met_p, RF_kern, u);
%fun = @(u) -sum(eval_likelihood_2D(img_p, met_p.alg, RF_kern, u), 'all');

% Maximize posterior probability for each frame
tic();
u_hat = fmincon(fun, u_0, [], [], [], [],...
       -0.5 * ones(size(u_0)),...
        0.5 * ones(size(u_0)), [], opt_p);

% Find best gain
a_gain = sum(u_hat .* u_mean, 'all') / norm(u_hat(:));

% Calculate error metrics (without gain)
err = error_metrics(u_hat, u_mean, 1);

fprintf(['i = %d | c_t =%5.2f | t =%4.0f\n' ...
    'err_0 = %.2f | err_1 = %.2f | A_opt = %.2f\n'], ...
    it, c_t, toc(), err_0(3), err(3), a_gain);

% Plot results
compare_frame(est_x, est_z, u_0 * 0.3, u_hat, u_mean)
