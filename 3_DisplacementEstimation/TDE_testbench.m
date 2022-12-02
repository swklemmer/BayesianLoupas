addpath('../lib/DispEst/')
addpath('../lib/SonoSim/BMS_aux/')

%% Testbench parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX');

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

% % Method parameters
met_param = struct(...
    'lambda',   1e3, ...    % prior weight
    'p',        1.05, ...      % norm deegree
    'vcn_z',    3, ...      % axial vecinity size [kernels]
    'vcn_x',    3 ...       % late. vecinity size [kernels]
    );         

% Material properties
ct_list = 0.25:0.25:3; % shear wave speed [m/s]
N_exp = 1;  % nr. of experiment per sws

for c_t = 1
    for n_i = 1:N_exp

        % Load results
        load(sprintf('../resources/BmodeData/ct%4.2f_%d.mat', ...
            c_t, n_i), 'RcvData', 'IData', 'QData', 'img_x', 'img_z')

        % Beamform IQ data [NOISE HAS'NT BEEN ADDED YET]
        [RF_mas, I_mas, Q_mas] = ...
            beamform(img_param, PData, Trans, TX, RcvData{1});

        % Split data into kernels
        [est_z, est_x, RF_kern, I_kern, Q_kern] = ...
            split_kernels(img_param, PData, RF_mas, I_mas, Q_mas);

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
graf = 2;

if graf == 1
    % Show Sonograms
    BFData = I_mas.^2 + Q_mas.^2;
    hObject.Value = 1;
    BMS_show_beamf(hObject, 0);
elseif graf == 2
    % Show estimations
    compare_estimates(est_x, est_z, u_0, u_hat)
end