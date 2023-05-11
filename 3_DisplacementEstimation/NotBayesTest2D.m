% This scripts estimates the displacements of a single US frame collection
% using a classical estimator (Loupas or NCC). Data is fetched from the
% "TrainFrames" folder.

addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')
addpath('../../poly2D/')
addpath(genpath('../../ImageWaves/'))

%% Simulation parameters

% Choose algorithm
alg = 'ncc';
oversample = 1; % Increase oversampling by a factor of d_samples [bool]
train_folder = 'SPW';

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_PWI.mat', ...
    'PData', 'Parameters', 'Trans', 'TX', 'TW', 'Receive');
lambda = 1540 / (TW(1).Parameters(1) * 1e6);

% Increase sample rate
if oversample
    train_folder = 'oversampled';
    d_samples = 20;
    PData.PDelta(3) = PData.PDelta(3) / d_samples;
    PData.Size(1) = d_samples * PData.Size(1);
end

% Imaging parameters
img_p = struct(...
    'f_c',      TW(1).Parameters(1) * 1e6, ...  % central frequency [Hz]
    'f_s',      TW(1).Parameters(1) / PData.PDelta(3) * 1e6, ...   % sampling frequency [Hz]
    't_s',      -1, ...                         % sampling frequency [Hz]
    'snr',      300 ...                         % signal-to-noise-ratio [dB]
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
    'crs_med',  0, ...      % coarse est. median filter
    'crs_sub',  .01, ...     % coarse sub-sample precission (ncc_poly only)
    'crs_win',  1, ...      % coarse est. windowing
    'fin_win',  0, ...      % fine est. windowing
    'moco_win', 0 ...       % motion correction windowing
   );

%% Perform estimation

% Load frames
c_t = 1.5; % [m/s]
it = 3;
load(sprintf('../resources/TrainData/%s/snr%d/ct%.2f_%d.mat', ...
    train_folder, img_p.snr, c_t, it), ...
    'RF_frames', 'I_frames', 'Q_frames', 'fr', 'c_t')

% Split data into kernels
[est_z, est_x, RF_kern, I_kern, Q_kern] = split_kernels(...
    est_p, PData, RF_frames, I_frames, Q_frames);

% Load true displacement
[u_tru, ~] = interp_real_u(est_x, est_z, fr, PData, lambda, c_t);

% Estimate displacemt using the chosen algorithm
tic()
if strcmp(alg, 'lou')
    u_hat = loupas_3D(est_p, I_kern, Q_kern);

elseif strcmp(alg, 'ncc')
    u_hat = ncc_3D(est_p, RF_kern) / img_p.f_s * img_p.f_c;

elseif strcmp(alg, 'ncc_poly')
    u_hat = ncc_poly_3D(est_p, RF_kern) / img_p.f_s * img_p.f_c; 
end

% Calculate error
err = error_metrics(u_hat, u_tru, 1);

% Find best gain
a_gain = sum(u_hat .* u_tru, 'all') / norm(u_hat(:));

fprintf([ ...
    'os. = %d | alg = %s\n' ...
    'i =%3d  | c_t =%5.2f | t =%4.1f\n' ...
    'G =%4.1f | e =%7.4f\n\n'], ...
    oversample, alg, it, c_t, toc(), err(3), a_gain);

% Plot estimation [wvls]
% figure()
% subplot(1, 2, 1)
% imagesc(est_x, est_z, u_hat, [-1, 1] * .1)
% subplot(1, 2, 2)
% imagesc(est_x, est_z, u_tru, [-1, 1] * .1)