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
    'z_len',    9, ...      % axial kernel length [wvls]
    'z_hop',    2, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    1, ...      % late. kernel hop    [wvls]
    't_len',    3 ...       % temp. kernel length [wvls]
    );

% % Method parameters
% met_param = struct(...
%     'lambda',   5e1, ...    % prior weight
%     'p',        1.05, ...   % norm deegree
%     'B',        5);         % vecinity size [kernels]


%%

% Load sonograms of different material properties
%ct_list = 0.25:0.25:3;  % shear wave speed [m/s]
ct_list = 1;
N_exp = 1;  % nr. of experiment per sws

for c_t = ct_list
    for n_i = 1:N_exp

        % Load results
        load(sprintf('../resources/BModeData/ct%4.2f_%d.mat', ...
            c_t, n_i), 'RcvData', 'IData', 'QData', 'img_x', 'img_z')

        % Beamform IQ data
        [RF_mas, I_mas, Q_mas] = beamform(img_param, PData, Trans, TX, RcvData{1});

        % MISSING: ADD NOISE

        % IF THE CURRENT METHOD DOESN'T WORK, WE CAN TRY TO ARTIFITIALLY
        % GENERATE RF AND IQ SIGNALS FROM THE FOLLOWING LINE
        %BFData = beamform_mag(img_param, IData{1}, QData{1});

        % Show Sonograms
%         BFData = RF_mas;
%         hObject.Value = 1;
%         BMS_show_beamf(hObject, 0);

        % Split data into kernels
        [RF_kern, I_kern, Q_kern] = ...
            split_kernels(img_param, PData, RF_mas, I_mas, Q_mas);

        % Calculate starting solution using Loupas
        u_0 = loupas_3D(I_kern, Q_kern);

%         % Maximize posterior probability
%         fun = @(u) -eval_posterior(img_param, met_param, rf_data, 'ack', u);
%         
%         opts = optimoptions('fmincon', ...
%             'Display', 'iter', ...
%             'MaxFunctionEvaluations', 5e4, ...
%             'OptimalityTolerance', 1e-4, ...
%             'StepTolerance', 1e-4);
%         u_hat = fmincon(fun, u_0, [], [], [], [],...
%             -0.5 * ones(img_param.K, 1),...
%             0.5 * ones(img_param.K, 1), [], opts);
    end
end