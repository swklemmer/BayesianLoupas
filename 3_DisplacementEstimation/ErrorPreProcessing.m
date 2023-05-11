% This script investigtes the effect of imaging system on displacement
% estimates and studies the preprocessing of the ground-truth-displacement
% data  before comparing it to estimated data.

addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')
addpath('../lib/TransCharac/')
addpath('../../Field_II/')
addpath(genpath('../../imagewaves/'))
addpath('../../../Vantage-4.7.6/Utilities/')

%% Simulation parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_PWI.mat', ...
    'PData', 'Parameters', 'Trans', 'TX', 'TW');
lambda = 1540 / (Trans.frequency * 1e6);

% Imaging parameters
img_p = struct(...
    'snr',      60, ...                     % signal-to-noise [dB]
    'f_c',      TW(1).Parameters(1) * 1e6, ...  % central frequency [Hz]
    'f_s',      40e6, ...                   % sampling frequency [Hz]
    't_s',      1 / 40e6, ...               % sampling frequency [Hz]
    'att',      0.3, ...                    % attenuation [dB/cm/MHz]
    'n_ang',    1, ...                      % steering angles
    'max_beta', 12, ...                     % max angle [º]
    'kaiser',   3);                         % Kaiser window beta

% Estimation parameters
est_p = struct(...
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    .5, ...     % axial kernel hop    [wvls]
    'x_len',    .5, ...     % late. kernel length [wvls]
    'x_hop',    .5, ...     % late. kernel hop    [wvls]
    't_len',    2, ...      % temp. kernel length [frames]
    'crs_mref', 1, ...      % coarse est. moving reference frame
    'fin_mref', 1, ...      % fine est. moving reference frame
    'crs_med',  1, ...      % coarse est. median filter
    'crs_win',  1, ...      % coarse est. windowing
    'moco_win', 1 ...       % motion correction windowing
   );

% Load estimated displacement data
ct_list = 0.25:0.25:3; % shear wave speed [m/s]
img_p.snr = 60; % [dB]

%% Estimate displacementsfrom many random RF realizations using Loupas

tic();
for c_t = 1:length(ct_list)
    for it = 1:10

        % Load a RF ensemble
        load(sprintf('../resources/TrainData/SPW/snr%d/ct%.2f_%d.mat', ...
            img_p.snr, ct_list(c_t), it), ...
            'I_frames', 'Q_frames', 'RF_frames', 'fr')
        
        % Split data into kernels
        [est_z, est_x, M_kern, I_kern, Q_kern] = ...
            split_kernels(est_p, PData, RF_frames, I_frames, Q_frames);
        
        % Pre-allocate estimated displacement
        if (c_t == 1) && (it == 1)
            u_est = zeros(length(est_z), length(est_x), length(ct_list), 10);
        end
        
        % Perform NCC estimation on it
        u_est(:, :, c_t, it) = loupas_3D(est_p, I_kern, Q_kern);

        fprintf('i=%2d | c_t=%4.2f | t=%2.0f\n', it, ct_list(c_t), toc());
    end
end

% Save estimations
mkdir(sprintf('../resources/EstData/lou/snr%d', img_p.snr));
save(sprintf('../resources/EstData/lou/snr%d/z%d.mat', ...
    img_p.snr, est_p.z_len), ...
    'u_est', 'est_p', 'est_x', 'est_z', 'ct_list', 'fr')

%% Load estimations, true displacement and obtain optimal correction param.

% Load estimations
load(sprintf('../resources/EstData/lou/snr%d/z%d.mat', ...
    img_p.snr, est_p.z_len), ...
    'u_est', 'est_p', 'est_x', 'est_z', 'ct_list', 'fr')

% Load true and peak displacement
u_tru = zeros(size(u_est, 1:3));
u_tru_peak = zeros(size(u_est, 3), 1);

for c_t = 1:length(ct_list)
    u_tru(:, :, c_t) = ...
        interp_real_u(est_x, est_z, fr, PData, lambda, ct_list(c_t));
    u_tru_peak(c_t) = max(u_tru(:, :, c_t), [], 'all');
end

% Calculate optimal gain for all estimations
a_gain = zeros(length(ct_list), 10);
a_lin = zeros(length(ct_list), 10);

for c_t = 1:length(ct_list)
    for it = 1:10
        a_gain(c_t, it) = ...
            sum(u_est(:, :, c_t, it) .* u_tru(:, :, c_t), 'all') / ...
            sum(u_est(:, :, c_t, it).^2, 'all');
    end
end

% Calculate error metrics before and after applying gain correction
u_err_0 = squeeze(rms(u_est - repmat(u_tru, [1 1 1 10]), [1 2 4]));
u_err_1 = squeeze(rms(...
    permute(repmat(a_gain, [1 1 size(u_est, [1 2])]), [3 4 1 2]) .* ...
    u_est - repmat(u_tru, [1 1 1 10]), [1 2 4]));

% Figures

fig = figure();
mesh(1:10, ct_list, a_gain)
xlabel('Realization')
ylabel('Shear wave speed [m/s]')
zlabel('Optimal gain')
zlim([0, 2])
view([110, 30])
title('Optimal gain for all realizations and SW speeds')
saveas(fig, sprintf('../results/img/optA1_z%d_snr%d', est_p.z_len, img_p.snr), 'jpg')

fig = figure();
plot(ct_list, mean(a_gain, 2))
xlabel('Shear wave speed [m/s]')
ylabel('Optimal gain')
ylim([0, 2])
title('Mean optimal gain for all SW speeds')
saveas(fig, sprintf('../results/img/optA2_z%d_snr%d', est_p.z_len, img_p.snr), 'jpg')

fig = figure();
plot(ct_list, 100 * u_err_0 ./ u_tru_peak)
hold on
plot(ct_list, 100 * u_err_1 ./ u_tru_peak)
hold off
ylabel('Peak %')
xlabel('c_t [m/s]')
legend({'Original', 'Gain corrected'})
ylim([0, 30])
title('Peak-normalized RMSE vs. shear wave speed')
saveas(fig, sprintf('../results/img/optA3_z%d_snr%d', est_p.z_len, img_p.snr), 'jpg')

%% Obtain local displacement deviation inside each res. cell
N_disp = 10;

% Setup libraries for resolution cell simulation
field_init(0);
Resource = struct('Parameters', Parameters);

% Obtain resolution cell dimensions [in wavelengths]
[~, x_wid_m, z_wid_m] = resolution_cell(img_p, Parameters, Trans, TW(1));
x_width = x_wid_m / lambda;
z_width = z_wid_m / lambda;
field_end();

% Load estimations
load(sprintf('../resources/EstData/lou/snr%d/z%d.mat', ...
    img_p.snr, est_p.z_len), ...
    'u_est', 'est_p', 'est_x', 'est_z', 'ct_list', 'fr')

% Load true displacement
u_tru = zeros(size(u_est, 1:3));

for c_t = 1:length(ct_list)
    u_tru(:, :, c_t) = interp_real_u( ...
        est_x, est_z, fr, PData, lambda, ct_list(c_t));
end

% Calculate error metrics
u_err_0 = squeeze(rms(u_est - repmat(u_tru, [1 1 1 10]), [1 2 4]));

% At each position, calculate deviation inside resolution cell
u_dev = zeros(length(est_z), length(est_x), length(ct_list));

tic();
for c_t = 1:length(ct_list)
    for z = 1:2:length(est_z)
        for x = 1:2:length(est_x)

            % Create new fine grid around current position
            fine_z = est_z(z) + z_width * ((0:N_disp-1) / N_disp - .5);
            fine_x = est_x(x) + x_width * ((0:N_disp-1) / N_disp - .5);

            % Obtain real displacement values at sampled points
            u_fine = interp_real_u(...
                fine_x, fine_z, fr, PData, lambda, ct_list(c_t));

            % Calculate their deviation
            u_dev(z, x, c_t) = std(u_fine, [], 'all');
        end
    end
    fprintf('c_t=%4.2f | t=%2.0f\n', ct_list(c_t), toc());
end

% Save results
save(sprintf('../results/stdev/snr%d_z%d.mat', img_p.snr, est_p.z_len), ...
    'N_disp', 'ct_list', 'est_p', 'img_p', 'dev_z', 'dev_x', 'u_dev', ...
    'x_width', 'x_wid_m', 'z_width', 'z_wid_m')


%% Correlate error and deviation as a function of location

z_len = 10;
snr = 60;

% Load real displacement deviation map
load(sprintf('../results/stdev/snr%d_z%d.mat', snr, z_len), ...
    'N_disp', 'est_p', 'img_p', 'dev_z', 'dev_x', 'u_dev', ...
    'x_width', 'x_wid_m', 'z_width', 'z_wid_m')

% Load estimations
load(sprintf('../resources/EstData/lou/snr%d/z%d.mat', ...
    img_p.snr, est_p.z_len), ...
    'u_est', 'est_p', 'est_x', 'est_z', 'ct_list', 'fr')

% Load true displacement
u_tru = interp_real_u(dev_x, dev_z, fr, PData, lambda, ct_list);

% Calculate error metrics
u_err_0 = abs(u_est(1:2:end, 1:2:end, :, :) - repmat(u_tru, [1 1 1 10]));


% Show scatter plot for each c_t 

for c_t = 1:length(ct_list)-1

    dev_pts = reshape(u_dev(:, :, c_t), [], 1);
    err_pts = reshape(u_err_0(:, :, c_t), [], 1);

    scatter(err_pts, dev_pts)
    pause()
end


%% Load ARF profile

% Filter real displacement using scaled ARF profile

% Calculate error metrics after preprocessing

% [Run optimization of filter width]

