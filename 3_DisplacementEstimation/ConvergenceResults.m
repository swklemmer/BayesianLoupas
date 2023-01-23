%% RMSE vs time

% Chose data to plot
ct = 1:4;
snr = 60;
metric = 3;
t_dim = 0:150;

% Load displacements
load(sprintf('../resources/TimeData/%s/alpha%.0f_%ddB.mat', ...
    'ack', 1, snr), 'err', 'elapsed_t')

% Interpolate times on fixed grid
ack1_e = zeros(length(ct), length(t_dim));

for j = ct

    % Find first zero in time vector
    last_t = find(elapsed_t(j, :) ~= 0, 1, 'last');

    % Last error is used as extrapolation value
    ack1_e(j, :) = ...
        interp1(elapsed_t(j, 1:last_t), err(j, 1:last_t, metric), ...
        t_dim, 'linear', err(j, last_t, metric));
end

% Average all frames
ack1_e = mean(ack1_e, 1);

load(sprintf('../resources/TimeData/%s/alpha%.0f_%ddB.mat', ...
    'ack', 2, snr), 'err', 'elapsed_t')

% Interpolate times on fixed grid
ack2_e = zeros(length(ct), length(t_dim));

for j = ct

    % Find first zero in time vector
    last_t = find(elapsed_t(j, :) ~= 0, 1, 'last');

    % Last error is used as extrapolation value
    ack2_e(j, :) = ...
        interp1(elapsed_t(j, 1:last_t), err(j, 1:last_t, metric), ...
        t_dim, 'linear', err(j, last_t, metric));
end

% Average all frames
ack2_e = mean(ack2_e, 1);

load(sprintf('../resources/TimeData/%s/alpha%.0f_%ddB.mat', ...
    'ncc', 1, snr), 'err', 'elapsed_t')

% Interpolate times on fixed grid
ncc1_e = zeros(length(ct), length(t_dim));

for j = ct

    % Find first zero in time vector
    last_t = find(elapsed_t(j, :) ~= 0, 1, 'last');

    % Last error is used as extrapolation value
    ncc1_e(j, :) = ...
        interp1(elapsed_t(j, 1:last_t), err(j, 1:last_t, metric), ...
        t_dim, 'linear', err(j, last_t, metric));
end

% Average all frames
ncc1_e = mean(ncc1_e, 1);

% Plot timeline
fig = figure(1);
fig.Position = [800, 200, 250, 200];
subplot(1, 1, 1)
plot(100 * t_dim / t_dim(end), ack1_e)
hold on
plot(100 * t_dim / t_dim(end), ack2_e)
plot(100 * t_dim / t_dim(end), ncc1_e)
hold off
grid on
ylim([6e-3, 11e-3])
xlim([1, 100])
xlabel('Optimization time [%]', 'FontSize', 12)
ylabel('RMSE [\lambda]', 'FontSize', 12)
sgtitle('Convergence time')
sgtitle(' ')
legend({'ACK L1', 'ACK L2', 'NCC L1'}, 'FontSize', 10, ...
    'Location', 'northeast')

% Convergence time
ack1_t = 1 + find(ack1_e(2:end) < (ack1_e(end) * 1.05), 1, 'first')
ack2_t = 1 + find(ack2_e(2:end) < (ack2_e(end) * 1.05), 1, 'first')
ncc1_t = 1 + find(ncc1_e(2:end) < (ncc1_e(end) * 1.05), 1, 'first')


%% Function evaluation rate

% Chose data to plot
snr = 60;

load(sprintf('../resources/TimeData/%s/alpha%.0f_%ddB.mat', ...
    'ack', 1, snr), 'fun_eval', 'elapsed_t')

ack1_eval_rate = mean(fun_eval(:, 1:end-1), 1) ./ mean(elapsed_t(:, 2:end), 1);

load(sprintf('../resources/TimeData/%s/alpha%.0f_%ddB.mat', ...
    'ack', 2, snr), 'fun_eval', 'elapsed_t')

ack2_eval_rate = mean(fun_eval(:, 1:end-1), 1) ./ mean(elapsed_t(:, 2:end), 1);

load(sprintf('../resources/TimeData/%s/alpha%.0f_%ddB.mat', ...
    'ncc', 1, snr), 'fun_eval', 'elapsed_t')

ncc1_eval_rate = mean(fun_eval(:, 1:end-1), 1) ./ mean(elapsed_t(:, 2:end), 1);


plot(ack1_eval_rate)
hold on
plot(ack2_eval_rate(1:7))
plot(ncc1_eval_rate)
hold off
legend({'ACK1', 'ACK 2', 'NCC 1'})

(mean(ack1_eval_rate) - mean(ncc1_eval_rate)) / mean(ncc1_eval_rate)
