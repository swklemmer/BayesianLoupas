addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')

%% Common parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_PWI.mat', ...
    'PData', 'Parameters', 'Trans', 'TX', 'TW', 'Receive');
lambda = 1540 / (TW(1).Parameters(1) * 1e6);
fr = 20;

%% Show estimations

alg = 'sqck';
p = 2;
snr = 60;
ct = 1;
ct_list = 1.25:0.5:3; % [m/s]
it = 2;

% Load estimations
load(sprintf('../resources/ErrorMetrics/%s-L%d/snr%ddB.mat', ...
    alg, p, snr), 'u_est')

% Load ground truth
load(sprintf('../resources/EstData/ground_truth/%.2f.mat', ...
    ct_list(ct)), 'u_mean', 'est_z', 'est_x', 'PData')

% Load fem displacement
[u_tru, ~] = interp_real_u(est_x, est_z, fr, PData, lambda, ct_list(ct));

% Plot displacements
iw(cat(3, permute(u_est(:, ct, it, :, :), [4 5 1 2 3]), 4* u_mean, u_tru), ...
    [-1, 1]*1e-1, 0)

%% Show estimation metrics

alg_list = {'sqck'};
p_list = [1 2];
snr_list = [60];
it_list = 1:4;

fig = figure(2);
fig.Position = [900, 0, 300, 700];
line_color = {'r', 'k', 'b', 'm'};
sgtitle('Effect of parameter \alpha on estimation error', 'FontSize', 14)
i = 0;

for alg = alg_list
    for p = p_list
        i = i + 1;
        for snr = snr_list
            % Load error metrics
            load(sprintf('../resources/ErrorMetrics/%s-L%d/snr%ddB.mat', ...
                alg{1}, p, snr), 'err', 'err_0', 'prm_list')

            % Obtain mean errors among all ct's and it's
            err_mean = squeeze(mean(err, 2:3));
            err0_mean = squeeze(mean(err_0(1, :, :, :), 2:3));

            % Bias
            subplot(3, 1, 1)
            semilogx(prm_list, err_mean(:, 1), line_color{i})
            yline(err0_mean(1), [line_color{i}, '--'])
            hold on
        
            % Std.
            subplot(3, 1, 2)
            semilogx(prm_list, sqrt(err_mean(:, 2)), line_color{i})
            yline(sqrt(err0_mean(2)), [line_color{i}, '--'])
            hold on
        
            % RMSE
            subplot(3, 1, 3)
            semilogx(prm_list, err_mean(:, 3), line_color{i})
            yline(err0_mean(3), [line_color{i}, '--'])
            hold on
        end
    end
end

ax1 = subplot(3, 1, 1);
hold off
grid on
ylabel('Bias [\lambda]', 'FontSize', 12)
%ylim([-4e-3 2e-3])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Location', 'northeast')

ax2 = subplot(3, 1, 2);
hold off
grid on
ylabel('St. Dev. [\lambda]', 'FontSize', 12)
%ylim([4e-3 12e-3])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Location', 'northwest')

ax3 = subplot(3, 1, 3);
hold off
grid on
ylabel('RMSE [\lambda]', 'FontSize', 12)
xlabel('Parameter \alpha', 'FontSize', 12)
%ylim([4e-3 12e-3])
legend({'5 dB', '', '20 dB', '', '60 dB', ''}, 'FontSize', 10, ...
    'Location', 'northwest')

linkaxes([ax1, ax2, ax3], 'x')

%% Compare algorithms
% THIS NEEDS FIXING
% THIS NEEDS FIXING
% THIS NEEDS FIXING

met_param.p = 1.05;
img_p.snr = 5;
ct_list = 2:4;
x_limits = [10^2.2, 10^5];

% Load Data
load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
'ack', met_param.p, img_p.snr), ...
'err', 'err0_mean', 'prm_list')

err_ack = squeeze(mean(err(:, ct_list, :), 2));
err0_mean = squeeze(mean(err0_mean(ct_list, :), 1));
x_ack = prm_list;

load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
'ncc', met_param.p, img_p.snr), ...
'err', 'prm_list')

err_ncc = squeeze(mean(err(:, ct_list, :), 2));
x_ncc = prm_list;

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
yline(err0_mean(1), [line_color{3}, '--'])
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
yline(sqrt(err0_mean(2)), [line_color{3}, '--'])
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
yline(err0_mean(3), [line_color{3}, '--'])
hold off
xlim(x_limits)
ylim([4e-3 12e-3])
ylabel('RMSE [\lambda]', 'FontSize', 12)
xlabel('Parameter \alpha', 'FontSize', 12)
legend({'ACK', 'NCC', 'Loupas'}, 'FontSize', 10, 'Location', 'southeast');

linkaxes([ax1, ax2, ax3], 'x')

%% Relative improvements

met_param.p = 1.05;
img_p.snr = 5;
ct_list = 1:4;

% Load Data
load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
    'ack', 2, img_p.snr), ...
    'err', 'err0_mean')

impr_ack = zeros(length(ct_list), 1);

for j = ct_list
    err0_mean = err0_mean(j, 3);
    err_ack = min(err(:, j, 3), [], 'all');

    impr_ack(j) = (err0_mean - err_ack) / err0_mean * 100;
end

impr_ack = mean(impr_ack)

load(sprintf('../resources/ErrorMetrics/%s/alpha%.0f_%ddB.mat', ...
    'ncc', met_param.p, img_p.snr), ...
    'err')

impr_ncc = zeros(length(ct_list), 1);

for j = ct_list
    err0_mean = err0_mean(j, 3);
    err_ncc = min(err(:, j, 3), [], 'all');

    impr_ncc(j) = (err0_mean - err_ncc) / err0_mean * 100;
end

%impr_ncc = mean(impr_ncc)