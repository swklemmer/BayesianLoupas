addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')

%% Imaging parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX', 'TW');
lambda = 1540 / (Trans.frequency * 1e6);

% Imaging parameters
img_param = struct(...
    'fr_n',     40, ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation (1 --> 20 kHz)
    'snr',      60, ...     % signal-to-noise-ratio [dB]
    'z_len',    8, ...      % axial kernel length [wvls]
    'z_hop',    2, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    2 ...       % temp. kernel length [wvls]
    );

% Method parameters
met_param = struct(...
    'alpha',    3, ...      % prior weight
    'p',        1.05, ...      % norm deegree
    'vcn_z',    3, ...      % axial vecinity size [kernels]
    'vcn_x',    1, ...      % late. vecinity size [kernels]
    'alg',      'ncc' ...   % Likelihood function
    );

%% Show displacement estimates

load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
    met_param.alg, met_param.p, img_param.snr), 'a_gain', 'u_est')

train_dir = sprintf('../resources/TrainData/SPW2/snr%d/', img_param.snr);
train_files = char(dir(fullfile(train_dir, '*.mat')).name);

for j = 1:size(train_files, 1)-3

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

    for i = 1:size(u_est, 1)

        % Load estimation
        u_hat = squeeze(u_est(i, j, :, :));
    
        % Show frames
        compare_frame(est_x, est_z, u_0, u_hat, u_tru)
        pause()
    end
end


%% Show estimation examples
img_param.snr = 60;
met_param.alg = 'ack';
met_param.p = 2;

load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
    met_param.alg, met_param.p, img_param.snr),...
    'err', 'err_0', 'param_list')

fig = figure(1);
fig.Position = [900, 200, 300, 300];
line_color = {'r', 'k', 'b', 'm'};
metric = 3;

for j = 2:size(err, 2)
    semilogx(param_list, err(:, j, metric), line_color{j})
    hold on
    yline(err_0(j, metric), [line_color{j}, '--'])
end

hold off
grid on
th = sgtitle('\alpha vs. estimation error', 'FontSize', 14);
ylabel('RMSE [\lambda]', 'FontSize', 12)
xlabel('Parameter \alpha', 'FontSize', 12)
ylim([3e-3 13e-3])
legend({'2 m/s', '', '2.5 m/s', '', '3 m/s', ''}, 'FontSize', 10, ...
    'Location', 'southwest')

%% Show estimation metrics

met_param.alg = 'ncc';
met_param.p = 1.05;
ct = 3;
snr_list = [5, 20, 60];

fig = figure(2);
fig.Position = [900, 800, 300, 700];
line_color = {'r', 'k', 'b', 'm'};
sgtitle('Effect of parameter \alpha on estimation error', 'FontSize', 14)

for i = 1:length(snr_list)
    load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
    met_param.alg, met_param.p, snr_list(i)), ...
    'err', 'err_0', 'param_list')
    
    err_mean = squeeze(mean(err(:, ct, :), 2));
    err_lou = squeeze(mean(err_0(ct, :), 1));

    % Bias
    subplot(3, 1, 1)
    semilogx(param_list, err_mean(:, 1), line_color{i})
    yline(err_lou(1), [line_color{i}, '--'])
    hold on

    % Std.
    subplot(3, 1, 2)
    semilogx(param_list, sqrt(err_mean(:, 2)), line_color{i})
    yline(sqrt(err_lou(2)), [line_color{i}, '--'])
    hold on

    % RMSE
    subplot(3, 1, 3)
    semilogx(param_list, err_mean(:, 3), line_color{i})
    yline(err_lou(3), [line_color{i}, '--'])
    hold on
end


ax1 = subplot(3, 1, 1);
hold off
grid on
ylabel('Bias [\lambda]', 'FontSize', 12)
ylim([-4e-3 2e-3])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Location', 'northeast')

ax2 = subplot(3, 1, 2);
hold off
grid on
ylabel('St. Dev. [\lambda]', 'FontSize', 12)
ylim([4e-3 12e-3])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Location', 'northwest')

ax3 = subplot(3, 1, 3);
hold off
grid on
ylabel('RMSE [\lambda]', 'FontSize', 12)
xlabel('Parameter \alpha', 'FontSize', 12)
ylim([4e-3 12e-3])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Location', 'northwest')

linkaxes([ax1, ax2, ax3], 'x')

%% Compare algorithms

met_param.p = 1.05;
img_param.snr = 5;
ct = 2:4;
x_limits = [10^2.2, 10^5];

% Load Data
load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
'ack', met_param.p, img_param.snr), ...
'err', 'err_0', 'param_list')

err_ack = squeeze(mean(err(:, ct, :), 2));
err_lou = squeeze(mean(err_0(ct, :), 1));
x_ack = param_list;

load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
'ncc', met_param.p, img_param.snr), ...
'err', 'param_list')

err_ncc = squeeze(mean(err(:, ct, :), 2));
x_ncc = param_list;

% Plot figure
fig = figure(3);
fig.Position = [900, 800, 300, 700];
line_color = {'r', 'k', 'b', 'm'};
sgtitle('Likelihood performance versus \alpha', 'FontSize', 14)

% Bias
ax1 = subplot(3, 1, 1);
semilogx(x_ack, err_ack(:, 1), line_color{1})
hold on
grid on
semilogx(x_ncc * 1e3, err_ncc(:, 1), line_color{2})
yline(err_lou(1), [line_color{3}, '--'])
hold off
xlim(x_limits)
ylim([-3e-3 2e-3])
ylabel('Bias [\lambda]', 'FontSize', 12)
legend({'ACK', 'NCC', 'Loupas'}, 'FontSize', 10, 'Location', 'northeast');

% Bias
ax2 = subplot(3, 1, 2);
semilogx(x_ack, sqrt(err_ack(:, 2)), line_color{1})
hold on
grid on
semilogx(x_ncc * 1e3, sqrt(err_ncc(:, 2)), line_color{2})
yline(sqrt(err_lou(2)), [line_color{3}, '--'])
hold off
xlim(x_limits)
ylim([4e-3 12e-3])
ylabel('St. Dev [\lambda]', 'FontSize', 12)
legend({'ACK', 'NCC', 'Loupas'}, 'FontSize', 10, 'Location', 'southeast');

% Bias
ax3 = subplot(3, 1, 3);
semilogx(x_ack, err_ack(:, 3), line_color{1})
hold on
grid on
semilogx(x_ncc * 1e3, err_ncc(:, 3), line_color{2})
yline(err_lou(3), [line_color{3}, '--'])
hold off
xlim(x_limits)
ylim([4e-3 12e-3])
ylabel('RMSE [\lambda]', 'FontSize', 12)
xlabel('Parameter \alpha', 'FontSize', 12)
legend({'ACK', 'NCC', 'Loupas'}, 'FontSize', 10, 'Location', 'southeast');

linkaxes([ax1, ax2, ax3], 'x')

%% Relative improvements

met_param.p = 1.05;
img_param.snr = 5;
ct = 1:4;

% Load Data
load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
    'ack', 2, img_param.snr), ...
    'err', 'err_0')

impr_ack = zeros(length(ct), 1);

for j = ct
    err_lou = err_0(j, 3);
    err_ack = min(err(:, j, 3), [], 'all');

    impr_ack(j) = (err_lou - err_ack) / err_lou * 100;
end

impr_ack = mean(impr_ack)

load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
    'ncc', met_param.p, img_param.snr), ...
    'err')

impr_ncc = zeros(length(ct), 1);

for j = ct
    err_lou = err_0(j, 3);
    err_ncc = min(err(:, j, 3), [], 'all');

    impr_ncc(j) = (err_lou - err_ncc) / err_lou * 100;
end

%impr_ncc = mean(impr_ncc)