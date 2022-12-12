SetUp_Sim;

%% Displacement estimation: 1PW vs 3PW DAS vs 3PW MAS Beamforming

% Simulate for different material properties
ct_list = 0.5:0.25:3;  % shear wave speed [m/s]
ct_list = 2.25;
N_exp = 1;              % nr. of experiment per sws
hObject.Value = 1;
P.bmode_adq = 38;

spw_error = zeros(length(ct_list) * N_exp, 3);
das_error = zeros(length(ct_list) * N_exp, 3);
mas_error = zeros(length(ct_list) * N_exp, 3);
mean_disp = zeros(length(ct_list) * N_exp, 1);
i = 1;

for c_t = ct_list
    for n_i = 1:N_exp
        tic();
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
            est_x, est_z, 1:size(MovieData, 3), PData, lambda, c_t);

        % Find optimal gains
        A_SPW = norm(SPW_est(:) - u_tru(:),  2) / norm(SPW_est(:),  2);
        A_DAS = norm(DAS_est(:) - u_tru(:),  2) / norm(DAS_est(:),  2);
        A_MAS = norm(MAS_est(:) - u_tru(:),  2) / norm(MAS_est(:),  2);

        % Calculate error metrics
        spw_error(i, :) = error_metrics(A_SPW * SPW_est, u_tru, mean_u);
        das_error(i, :) = error_metrics(A_DAS * DAS_est, u_tru, mean_u);
        mas_error(i, :) = error_metrics(A_MAS * MAS_est, u_tru, mean_u);
        mean_disp(i) = mean_u;
        i = i + 1;

        fprintf('c_t = %4.2f | n = %d | t = %3.0f\n' , c_t, n_i, toc());

        % Show displacements
        compare_estimates(est_x, est_z, ...
            A_SPW * SPW_est, A_DAS * DAS_est, u_tru)
    end
end

% Process Errors
spw_mean = mean(spw_error(11:end, :), 1);
spw_std = std(spw_error(11:end, :), [], 1);
das_mean = mean(das_error(11:end, :), 1);
das_std = std(das_error(11:end, :), [], 1);
mas_mean = mean(mas_error(11:end, :), 1);
mas_std = std(mas_error(11:end, :), [], 1);

% Save
save('../resources/ErrorMetrics/ImgSetup.mat', 'ct_list', 'N_exp', ...
    'spw_error', 'das_error', 'mas_error', 'mean_disp')
