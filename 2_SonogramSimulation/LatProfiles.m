SetUp_Sim;

%% Lateral profiles
c_t = 2;
n_i = 1;

% Load sonograms
load(sprintf('../resources/BModeData/SPW/ct%4.2f_%d.mat', ...
    c_t, n_i), 'img_x', 'img_z', 'IData', 'QData')

% Estimate displacement using Single Plane Wave
param_flag = 1;
BMS_estimate_u_lou(IData{1}, QData{1});
SPW_est = MovieData;

% Load sonograms
load(sprintf('../resources/BModeData/MAS/ct%4.2f_%d.mat', ...
    c_t, n_i), 'RcvData')

% Beamform IQ data using L-DAS beamforming
[~, I_das, Q_das] = ...
    beamform_rf_lin(img_param, PData, Trans, TX, RcvData{1});

% Estimate displacement using L-DAS beamforming
I_das = reshape(I_das, [size(I_das, 1, 2), 1, size(I_das, 3)]);
Q_das = reshape(Q_das, [size(Q_das, 1, 2), 1, size(Q_das, 3)]);
param_flag = 1;
BMS_estimate_u_lou(I_das, Q_das);
DAS_est = MovieData;

% Beamform IQ data using NL-MAS beamforming
[~, I_mas, Q_mas] = ...
    beamform_rf(img_param, PData, Trans, TX, RcvData{1});

% Estimate displacement using NL-MAS beamforming
I_mas = reshape(I_mas, [size(I_mas, 1, 2), 1, size(I_mas, 3)]);
Q_mas = reshape(Q_mas, [size(Q_mas, 1, 2), 1, size(Q_mas, 3)]);
param_flag = 1;
BMS_estimate_u_lou(I_mas, Q_mas);
MAS_est = MovieData;

% Load true displacement
[u_tru, mean_u]= interp_real_u(...
    est_x, est_z, size(MovieData, 3), lambda, PData, c_t);

% Find optimal gains
A_SPW = norm(SPW_est(:) - u_tru(:),  2) / norm(SPW_est(:),  2);
A_DAS = norm(DAS_est(:) - u_tru(:),  2) / norm(DAS_est(:),  2);
A_MAS = norm(MAS_est(:) - u_tru(:),  2) / norm(MAS_est(:),  2);

% Plot lateral profiles

z_plot = 22;
%t_plot = 13;
for t_plot = 5:36
    plot(est_x, A_SPW * SPW_est(z_plot, :, t_plot))
    hold on
    plot(est_x, A_DAS * DAS_est(z_plot, :, t_plot), 'o-')
    plot(est_x, A_MAS * MAS_est(z_plot, :, t_plot), '*')
    plot(est_x, u_tru(z_plot, :, t_plot - 4), '--')
    plot(est_x, u_tru(z_plot-1, :, t_plot - 4), '--')
    plot(est_x, u_tru(z_plot+1, :, t_plot - 4), '--')
    hold off
    grid on
    ylim([-0.1, 0.3])
    legend({'SPW', 'DAS', 'MAS', 'true'})
    xlabel('Lateral dist. [\lambda]')
    ylabel('Displacement [\lambda]')
    pause()
end
