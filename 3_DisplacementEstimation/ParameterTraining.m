addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')
graf = 0;

%% Training parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_PWI.mat', ...
    'PData', 'Parameters', 'Trans', 'TX', 'TW', 'Receive');
lambda = 1540 / (TW(1).Parameters(1) * 1e6);

% Imaging parameters
img_p = struct(...
    'f_c',      TW(1).Parameters(1) * 1e6, ...  % central frequency [Hz]
    'f_s',      TW(1).Parameters(1) / PData.PDelta(3) * 1e6, ... % s.f.[Hz]
    't_s',      -1, ...                         % sampling frequency [Hz]
    'snr',      -1 ...                         % signal-to-noise-ratio [dB]
    );

img_p.t_s = 1 / img_p.f_s;

% Estimation parameters
est_p = struct(...
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    2, ...     % axial kernel hop    [wvls]
    'x_len',    .5, ...    % late. kernel length [wvls]
    'x_hop',    2, ...     % late. kernel hop    [wvls]
    't_len',    2, ...     % temp. kernel length [frames]
    'fin_mref', 1, ...     % fine est. moving reference frame
    'fin_win',  1 ...      % fine est. windowing
    );

% Method parameters
met_p = struct(...
    'alg',     -1, ...  % Likelihood function
    'alpha',   -1, ...  % prior weight
    'p',       -1, ...  % norm deegree
    'vcn_z',    3, ...  % axial vecinity size [kernels]
    'vcn_x',    1 ...   % late. vecinity size [kernels]
    );

% Optimization parameters
opt_p = optimoptions('fmincon', ...
    'Display', 'none', ...
    'MaxFunctionEvaluations', 2e4, ...
    'MaxIterations', 10, ...
    'OptimalityTolerance', 1e-2, ...
    'StepTolerance', 1e-5 ...
    );

%% Calculate errors
alg_list = {'ack', 'ncc', 'sqck'};
p_list = [1, 2];
snr_list = [5 20 60];
prm_list = logspace(0, 2, 10);
ct_list = 1.25:0.5:3; % [m/s]
it_list = 5:7;


% Pre-allocate error metrics
err = zeros(length(prm_list), length(ct_list), length(it_list), 3);
err_0 = zeros(length(prm_list), length(ct_list), length(it_list), 3);
elap_t = zeros(length(prm_list), length(ct_list), length(it_list));

tic()
for alg = alg_list
    for p = p_list
        met_p.p = p;
        met_p.alg = alg{1};

        % Choose alpha range
        if (p == 1)
            if strcmp(met_p.alg, 'ack')
                prm_list = logspace(0, 2, 10);
            elseif strcmp(met_p.alg, 'ncc')
                prm_list = logspace(-1.5, 0.5, 10);
            elseif strcmp(met_p.alg, 'sqck')
                prm_list = logspace(0, 2, 10);
            end
        elseif (p == 2)
            if strcmp(met_p.alg, 'ack')
                prm_list = logspace(-0.5, 1, 10);
            elseif strcmp(met_p.alg, 'ncc')
                prm_list = logspace(-2, -1, 10);
            elseif strcmp(met_p.alg, 'sqck')
                prm_list = logspace(-2, 0, 10);
            end
        end

        for snr_i = 1:length(snr_list)
            img_p.snr = snr_list(snr_i);

            for prm_j = 1:length(prm_list)

                % Raise parameter change flag
                param_flag = 1;

                % Update parameter
                met_p.alpha = prm_list(prm_j);

                for ct_k = 1:length(ct_list)

                    % Load ground truth
                    load(sprintf( ...
                        '../resources/EstData/ground_truth/%.2f.mat', ct_list(ct_k)), ...
                        'u_mean');

                    for it_l = it_list
                        % Load training data
                        load(sprintf('../resources/TrainData/SPW/snr%d/ct%.2f_%d', ...
                            snr_list(snr_i), ct_list(ct_k), it_l), ...
                            'RF_frames', 'I_frames', 'Q_frames', 'fr')

                        % Split data into kernels
                        [est_z, est_x, RF_kern, I_kern, Q_kern] = ...
                            split_kernels(est_p, PData, RF_frames, I_frames, Q_frames);

                        % Calculate starting solution using Loupas and its error metrics
                        u_0 = loupas_3D(est_p, I_kern, Q_kern);
                        err_0(prm_j, ct_k, it_l, :) = error_metrics(u_0, u_mean, 1);

                        % Pre-compute maximum fast freqeuncy
                        max_fc = pre_comp_SQCK(img_p, u_0, RF_kern);

                        if (prm_j == 1) && (ct_k == 1) && (it_l == 1)
                            % Pre-allocate estimations
                            u_est = zeros([length(prm_list), length(ct_list), length(it_list), size(u_0)]);
                        end

                        % Declare loss function
                        fun = @(u) -eval_posterior_2D(img_p, met_p, RF_kern, u);

                        % Maximize posterior probability for each frame
                        u_hat = fmincon(fun, u_0, [], [], [], [],...
                            -0.5 * ones(size(u_0)),...
                            0.5 * ones(size(u_0)), [], opt_p);

                        % Calculate error metrics (without gain)
                        err(prm_j, ct_k, it_l, :) = error_metrics(u_hat, u_mean, 1);

                        fprintf('alg = %s-L%d | param = %d | c_t =%5.2f | it = %d | t =%4.0f\n', ...
                            alg{1}, p, prm_j, ct_list(ct_k), it_l, toc());

                        % Save estimation
                        u_est(prm_j, ct_k, it_l, :, :) = u_hat;

                    end
                end
            end

            if ~graf
                mkdir(sprintf('../resources/ErrorMetrics/comp_%s-L%d', met_p.alg, met_p.p))
                save(sprintf('../resources/ErrorMetrics/comp_%s-L%d/snr%ddB.mat', ...
                    met_p.alg, met_p.p, img_p.snr), ...
                    'est_x', 'est_z', 'u_est', 'err', 'err_0', ...
                    'prm_list', 'img_p', 'met_p', 'opt_p')
            end
        end
    end

end
