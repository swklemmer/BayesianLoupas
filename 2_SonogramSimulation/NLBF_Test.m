addpath('../lib/SonoSim/');
addpath('../lib/SonoSim/BMS_aux/');
addpath('../lib/DispEst/')

%% Load simulation parameters

% Specify user defined parameters
P = struct(...
    'startDepth',   25, ...     % acq. start depth [wvls]
    'endDepth',     75, ...     % acq. end depth [wvls]
    'latDist',      30, ...     % acq. lateral span [wvls]
    'bmode_dly',    50, ...    % delay between b-mode images [usec] (20kHz)
    'bmode_adq',    2e3 / 50, ... % b-mode acq. number
    'hv',           1.6, ...    % transmition bipolar voltage
    'n_ang',        3, ...      % nr. of steering angles
    'n_push',       32, ...     % nr. of pushing elements
    'z_push',       10e-3, ...  % Push depth [m]
    'c_push',       430, ...    % nr. of push cycles
    'simulate',     1, ...      % enable simulate mode
    'c_t',          1.5 ...     % simulated shear wave speed [m/s]
    );

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX');
lambda = 1540 / (Trans.frequency * 1e6);

% Imaging parameters
img_param = struct(...
    'fr_n',     40, ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation (1 --> 20 kHz)
    'snr',      20, ...     % signal-to-noise-ratio [dB]
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    5, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    3 ...       % temp. kernel length [wvls]
    );

% Estimation parameters
param_flag = 1;
current_param = struct( ...
    'ncc', struct( ...
        'axi_len', 9, ...       % NCC: Axial window length [wvls]
        'axi_hop', 1, ...       % NCC: Axial window hop [wvls]
        'lat_len', 1, ...       % NCC: Lateral window length [wvls]
        'lat_hop', 1, ...       % NCC: Lateral window hop [wvls]
        'fine_res', 0.1, ...    % NCC: Polynomial interp. res. [smpls]
        'med_sz', 7), ...       % NCC: Median filter size [smpls]
    'scc', struct( ...
        'axi_len', 9, ...       % SCC: Axial window length [wvls]
        'axi_hop', 1, ...       % SCC: Axial window hop [wvls]
        'lat_len', 1, ...       % SCC: Lateral window length [wvls]
        'lat_hop', 1, ...       % SCC: Lateral window hop [wvls]
        'fine_res', 0.1, ...    % SCC: Polynomial interp. res. [smpls]
        'med_sz', 5, ...        % SCC: Median filter size [smpls]
        'search_z', 2, ...      % SCC: Axial disp. limit [wvls]
        'search_x', 2), ...     % SCC: Lateral disp. limit [wvls]
    'lou', struct( ...
        'axi_len', 9, ...       % LOU: Axial window length [wvls]
        'axi_hop', 1, ...       % LOU: Axial window hop [wvls]
        'lat_len', 1, ...       % LOU: Lateral window length [wvls]
        'lat_hop', 1, ...       % LOU: Lateral window hop [wvls]
        'med_sz', 7, ...        % LOU: Median filter size [smpls]
        'ens_len', 3, ...       % LOU: N = Ensemble length
        'cum_sum', 15));    % LOU: Moving average size [smpls]

% Simulate for different material properties
ct_list = 0.25:0.75:3;  % shear wave speed [m/s]
N_exp = 1;              % nr. of experiment per sws
hObject.Value = 1;

%% Displacement estimation: DAS vs MAS Beamforming
P.bmode_adq = 38;

spw_error = zeros(length(ct_list), 3);
das_error = zeros(length(ct_list), 3);
mas_error = zeros(length(ct_list), 3);
max_disp = zeros(length(ct_list), 1);

for c = 1:length(ct_list)
    for n_i = 1:N_exp
        % Load sonograms
        load(sprintf('../resources/BModeData/SPW/ct%4.2f_%d.mat', ...
            ct_list(c), n_i), 'img_x', 'img_z', 'IData', 'QData')

        % Estimate displacement using Single Plane Wave
        param_flag = 1;
        BMS_estimate_u_lou(IData{1}, QData{1});
        SPW_est = MovieData;

        % Load sonograms
        load(sprintf('../resources/BModeData/MAS/ct%4.2f_%d.mat', ...
            ct_list(c), n_i), 'RcvData')

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
        [u_tru, max_u]= interp_real_u(...
            est_x, est_z, size(MovieData, 3), lambda, PData, ct_list(c));

        % Find optimal gains
        A_SPW = norm(SPW_est(:) - u_tru(:),  2) / norm(SPW_est(:),  2);
        A_DAS = norm(DAS_est(:) - u_tru(:),  2) / norm(DAS_est(:),  2);
        A_MAS = norm(MAS_est(:) - u_tru(:),  2) / norm(MAS_est(:),  2);

        % Calculate error metrics
        spw_error(c, :) = error_metrics(A_SPW * SPW_est, u_tru);
        das_error(c, :) = error_metrics(A_DAS * DAS_est, u_tru);
        mas_error(c, :) = error_metrics(A_MAS * MAS_est, u_tru);
        max_disp(c) = max_u;

        % Show displacements
%         compare_estimates(est_x, est_z, ...
%             A_DAS * DAS_est, A_MAS * MAS_est, u_tru)
    end
end

%% NL-Beamforming over Magnitude Data

for c = flip(ct_list)    
    for n_i = 1:N_exp
        % Load results
        load(sprintf('../resources/BModeData/MAS/ct%4.2f_%d.mat', ...
            c, n_i), 'img_x', 'img_z', 'BFData')

        % Show displacement
        BMS_show_beamf(hObject, 0, c);
    end
end

%% NL-Beamforming over RF Data

for c = flip(ct_list)    
    for n_i = 1:N_exp
        % Load results
        load(sprintf('../resources/BmodeData/MAS/ct%4.2f_%d.mat', ...
            c, n_i), 'img_x', 'img_z', 'RcvData')

        % Beamform IQ data
        [RF_mas, I_mas, Q_mas] = ...
            beamform_rf_lin(img_param, PData, Trans, TX, RcvData{1});

        % Show displacement
        BFData = I_mas.^2 + Q_mas.^2;
        BMS_show_beamf(hObject, 0, c);
    end
end

%% NL-Beamforming over IQ Data

for c = flip(ct_list)    
    for n_i = 1:N_exp
        % Load results
        load(sprintf('../resources/BModeData/SPW/ct%4.2f_%d.mat', ...
            c, n_i), 'img_x', 'img_z', 'IData', 'QData')

        % Show displacement
        BFData = squeeze(IData{1}.^2 + QData{1}.^2);
        BMS_show_beamf(hObject, 0, c);
    end
end