addpath('../lib/SonoSim/')
addpath('../lib/DispEst/')
graf = 0;

%% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX', 'TW');

% Imaging parameters
img_param = struct(...
    'fr_n',     40, ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation (1 --> 20 kHz)
    'snr',      60, ...     % signal-to-noise-ratio [dB]
    'n_ang',    1, ...      % steering angles
    't_len',    2 ...       % temp. kernel length [wvls]
    );

snr_list = [5:5:30, 60]; % signal-to-noise-ratio [dB]
ct_list = 1.5:0.5:3; % shear wave speed [m/s]

% Save directory
save_dir = 'SPW2';

RF_frames = zeros([PData.Size([1 2]), 3]);
I_frames = zeros([PData.Size([1 2]), 3]);
Q_frames = zeros([PData.Size([1 2]), 3]);
rng(6942)

for snr = snr_list
    img_param.snr = snr;
    mkdir(sprintf('../resources/TrainData/%s/snr%d', save_dir, img_param.snr))
    for c_t = ct_list 
        tic();
        % Load bmode images
        load(sprintf('../resources/BmodeData/SPW/ct%4.2f_1.mat', c_t), ...
            'RcvData', 'img_x', 'img_z')

        % Select center frame at random
        fr = randi([15, img_param.fr_n - img_param.t_len]);
        frames = fr + (1:img_param.t_len) - ceil(img_param.t_len/2);

        % Beamform RF data at selected frames
        [RF_frames, I_frames, Q_frames] = ...
            beamform_rf_lin(img_param, PData, Trans, TX, RcvData{1}, frames);
    
        % Save beamformed frame
        save(sprintf('../resources/TrainData/%s/snr%d/ct%4.2f.mat', ...
        save_dir, img_param.snr, c_t), ...
        'fr', 'c_t', 'RF_frames', 'I_frames', 'Q_frames')

        fprintf('%ddB | c_t = %4.2f | t = %3.0f\n', img_param.snr, c_t, toc());

        % Show center B-mode frame
        if graf
            figure(1)
            Mag_frame = I_frames(:, :, floor(end/2)).^2 + ...
                        Q_frames(:, :, floor(end/2)).^2;
            imagesc(img_x, img_z, Mag_frame)
            colorbar
        end
    end
end