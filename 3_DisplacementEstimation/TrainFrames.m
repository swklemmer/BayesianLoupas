addpath('../lib/SonoSim/')
addpath('../lib/DispEst/')
graf = 0;

%% Simulation parameters

oversample = 1; % Increase oversampling by a factor of d_samples [bool]

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_PWI.mat', ...
    'PData', 'Parameters', 'Trans', 'TX', 'TW', 'Receive');
lambda = 1540 / (TW(1).Parameters(1) * 1e6);

% Increase sample rate
if oversample
    d_samples = 20;
    PData.PDelta(3) = PData.PDelta(3) / d_samples;
    PData.Size(1) = d_samples * PData.Size(1);
end

% Imaging parameters
img_p = struct(...
    'f_c',      TW(1).Parameters(1) * 1e6, ...  % central frequency [Hz]
    'f_s',      TW(1).Parameters(1) / PData.PDelta(3) * 1e6, ... % sam. fr
    't_s',      -1, ...         % sampling frequency [Hz]
    'snr',      300, ...        % signal-to-noise-ratio [dB]
    'n_ang',    1, ...          % nr. of steering angles
    'fr_n',     40 ...          % nr. of input frames
    );

img_p.t_s = 1 / img_p.f_s;

% Estimation parameters
est_p = struct(...
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    2, ...      % axial kernel hop    [wvls]
    'x_len',    .5, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    2, ...      % temp. kernel length [frames]
    'crs_mref', 0, ...      % coarse est. moving reference frame
    'fin_mref', 1, ...      % fine est. moving reference frame
    'crs_med',  1, ...      % coarse est. median filter
    'crs_sub',  .1, ...     % coarse sub-sample precission (ncc_poly only)
    'crs_win',  1, ...      % coarse est. windowing
    'fin_win',  0, ...      % fine est. windowing
    'moco_win', 0 ...       % motion correction windowing
   );

snr_list = [60 20 5]; % signal-to-noise-ratio [dB]
ct_list = 1.25:0.5:3; % shear wave speed [m/s]

% Save directory
save_dir = '../resources/TrainData/SPW';
mkdir(save_dir);

%% Start simulation

RF_frames = zeros([PData.Size([1 2]), 3]);
I_frames = zeros([PData.Size([1 2]), 3]);
Q_frames = zeros([PData.Size([1 2]), 3]);

tic();
for snr = snr_list
    img_p.snr = snr;
    mkdir(sprintf('%s/snr%d', save_dir, img_p.snr))
    for c_t = ct_list 
        for it = 1:4
            % Load bmode images
            load(sprintf('../resources/BmodeData/SPW/ct%4.2f_%d.mat', ...
                c_t, it), 'RcvData', 'img_x', 'img_z')
    
            % Select center frame
            fr = 20;
            frames = fr + (1:est_p.t_len) - ceil(est_p.t_len/2);
    
            % Beamform RF data at selected frames
            [RF_frames, I_frames, Q_frames] = beamform_rf_lin(...
                img_p, PData, Trans, TX, Receive, RcvData{1}, frames);
        
            % Save beamformed frame
            save(sprintf('%s/snr%d/ct%4.2f_%d.mat', ...
                save_dir, snr, c_t, it), ...
                'fr', 'c_t', 'RF_frames', 'I_frames', 'Q_frames')
    
            fprintf('%ddB | c_t = %4.2f | it = %d | t = %3.0f\n', ...
                img_p.snr, c_t, it, toc());
    
            % Show center B-mode frame
            if graf
                figure(1)
                Mag_frame = I_frames(:, :, floor(end/2)).^2 + ...
                            Q_frames(:, :, floor(end/2)).^2;
                imagesc(Mag_frame)
                colorbar
            end
        end
    end
end