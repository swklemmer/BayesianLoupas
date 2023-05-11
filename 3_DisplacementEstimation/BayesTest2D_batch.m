% This scripts estimates the displacements of a batch of US frames
% corresponding to the full propagation of a shear wave using bayesian
% estimation (may take very long if input frame number is high!). 

addpath('../lib/DispEst/')
addpath('../lib/SonoSim/BMS_aux/')

%% Testbench parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_PWI.mat', ...
    'PData', 'Parameters', 'Trans', 'TX', 'TW', 'Receive');
lambda = 1540 / (TW(1).Parameters(1) * 1e6);

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
    'z_hop',    1, ...     % axial kernel hop    [wvls]
    'x_len',    .5, ...    % late. kernel length [wvls]
    'x_hop',    .5, ...     % late. kernel hop    [wvls]
    't_len',    2, ...     % temp. kernel length [frames]
    'fin_mref', 1, ...     % fine est. moving reference frame
    'fin_win',  1 ...      % fine est. windowing
   );

% Method parameters
met_param = struct(...
    'lambda',   1e3, ...    % prior weight
    'p',        1.05, ...      % norm deegree
    'vcn_z',    3, ...      % axial vecinity size [kernels]
    'vcn_x',    3 ...       % late. vecinity size [kernels]
    );         

% Material properties
ct_list = 0.75:0.25:3; % shear wave speed [m/s]
N_exp = 1:10;  % nr. of experiment per sws

for c_t = 1
    for n_i = 1:N_exp

        % Load results
        load(sprintf('../resources/BmodeData/ct%4.2f_%d.mat', ...
            c_t, n_i), 'RcvData', 'IData', 'QData', 'img_x', 'img_z')

        % Beamform RF data at selected frames
        [RF_frames, I_frames, Q_frames] = beamform_rf_lin(...
            img_p, PData, Trans, TX, Receive, RcvData{1});

        % Split data into kernels
        [est_z, est_x, RF_kern, I_kern, Q_kern] = ...
            split_kernels(img_param, PData, RF_frames, I_frames, Q_frames);

        % Calculate starting solution using Loupas
        u_0 = loupas_3D(I_kern, Q_kern);

        % Maximize posterior probability for each frame
        opts = optimoptions('fmincon', ...
            'Display', 'iter-detailed', ...
            'MaxFunctionEvaluations', 2e4, ...
            'MaxIterations', 3e1, ...
            'OptimalityTolerance', 1e-4, ...
            'StepTolerance', 1e-4);

        u_hat = zeros(size(u_0));

        for t = 1:size(u_0, 3)
            % Select kernel corresponding to current frame
            fun = @(u) -eval_posterior_2D(met_param,...
                RF_kern(:, :, t, :, :, :), 'ack', u);

            u_hat(:, :, t) = fmincon(fun, u_0(:, :, t), [], [], [], [],...
               -0.5 * ones(size(u_0, 1, 2)),...
                0.5 * ones(size(u_0, 1, 2)), [], opts);
        end
    end
end


%% Animations
graf = 1;

if graf == 1
    % Show Sonograms
    BFData = I_frames.^2 + Q_frames.^2;
    hObject.Value = 1;
    BMS_show_beamf(hObject, 0);
elseif graf == 2
    % Show estimations
    compare_estimates(est_x, est_z, u_0, u_hat)
end