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
    'snr',      5, ...     % signal-to-noise-ratio [dB]
    'z_len',    8, ...      % axial kernel length [wvls]
    'z_hop',    2, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    2 ...       % temp. kernel length [wvls]
    );

%% Show displacement estimates

load(sprintf('../resources/ErrorMetrics/alpha2_%ddB.mat', img_param.snr), ...
'a_gain', 'u_est')

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

%         % Calculate error metrics (without gain)
%         err(i, j, :) = error_metrics(u_hat, u_tru, 1);

        pause()
    end
end


%% Show estimation examples
img_param.snr = 60;

load(sprintf('../resources/ErrorMetrics/alpha1_%ddB.mat', img_param.snr), ...
'err', 'err_0', 'param_list')

fig = figure(1);
fig.Position = [900, 200, 500, 300];
line_color = {'r', 'k', 'b', 'm'};
metric = 3;

for j = 2:size(err, 2)
    semilogx(param_list, err(:, j, metric), line_color{j})
    hold on
    yline(err_0(j, metric), [line_color{j}, '--'])
end

hold off
grid on
title('Effect of parameter \alpha on estimation error.', 'FontSize', 14)
ylabel('RMSE [\lambda]', 'FontSize', 12)
xlabel('Parameter \alpha', 'FontSize', 12)
ylim([2.5e-3 12.5e-3])
legend({'2 m/s', '', '2.5 m/s', '', '3 m/s', ''}, 'FontSize', 12, 'Position', [0.17 0.16 0.1 0.1])

%% Show estimation metrics
snr_list = [5, 20, 60];

fig = figure(2);
fig.Position = [900, 800, 400, 700];
line_color = {'r', 'k', 'b', 'm'};

for i = 1:length(snr_list)
    load(sprintf('../resources/ErrorMetrics/alpha2_%ddB.mat', snr_list(i)), ...
    'err', 'err_0', 'param_list')
    
    err_mean = squeeze(mean(err(:, [1:4], :), 2));

    % Bias
    subplot(3, 1, 1)
    semilogx(param_list, err_mean(:, 1), line_color{i})
    yline(err_0(3, 1), [line_color{i}, '--'])
    hold on

    % Std.
    subplot(3, 1, 2)
    semilogx(param_list, sqrt(err_mean(:, 2)), line_color{i})
    yline(sqrt(err_0(3, 2)), [line_color{i}, '--'])
    hold on

    % RMSE
    subplot(3, 1, 3)
    semilogx(param_list, err_mean(:, 3), line_color{i})
    yline(err_0(3, 3), [line_color{i}, '--'])
    hold on
end

sgtitle('Effect of parameter \alpha on estimation error.', 'FontSize', 14)

subplot(3, 1, 1)
hold off
grid on
ylabel('Bias [\lambda]', 'FontSize', 12)
ylim([-4e-3 1e-3])
%xlim([0.25 20])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Position', [0.17 0.7 0.1 0.05])

subplot(3, 1, 2)
hold off
grid on
ylabel('St. Dev. [\lambda]', 'FontSize', 12)
ylim([4e-3 12e-3])
%xlim([0.25 20])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Position', [0.17 0.407 0.1 0.05])

subplot(3, 1, 3)
hold off
grid on
ylabel('RMSE [\lambda]', 'FontSize', 12)
xlabel('Parameter \alpha', 'FontSize', 12)
ylim([4e-3 12e-3])
%xlim([0.25 20])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Position', [0.17 0.115 0.1 0.05])
