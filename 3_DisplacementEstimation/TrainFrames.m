addpath('../lib/SonoSim/')
addpath('../lib/DispEst/')
graf = 0;

%% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX', 'TW');

% Imaging parameters
img_param = struct(...
    'fr_n',     40, ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation (1 --> 20 kHz)
    'snr',      5, ...       % signal-to-noise-ratio [dB]
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    5, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    3 ...       % temp. kernel length [wvls]
    );

snr_list = [5:5:30, 60]; % signal-to-noise-ratio [dB]
ct_list = 0.25:0.75:3; % shear wave speed [m/s]

RF_frames = zeros([PData.Size([1 2]), 3]);
I_frames = zeros([PData.Size([1 2]), 3]);
Q_frames = zeros([PData.Size([1 2]), 3]);
rng(69420)

for snr = snr_list
    snr
    img_param.snr = snr;
    for c_t = ct_list
        c_t
    
        % Load bmode images
        load(sprintf('../resources/BmodeData/ct%4.2f_1.mat', c_t), ...
            'BFData', 'RcvData', 'IData', 'QData', 'img_x', 'img_z')
    
        % Beamform RF data
        [RF_mas, I_mas, Q_mas] = ...
            beamform_rf(img_param, PData, Trans, TX, RcvData{1});
    
        % Load displacement info from .h5 file
        [rea_t, rea_x, rea_y, rea_z, ~, ~, u_z] = ...
            load_u(sprintf('../resources/FemData/u_%.0f.h5', 3e3 * c_t^2));
        
        % Frame to time transformation
        fr2t = 50 * 1e-6 / diff(rea_t([1, 2]));
    
        % Pick a center frame at random
        fr = 1 + randi([2, size(RF_mas, 3) - 1]);
    
        % Obtain 3 frames from beamformed data
        RF_frames = RF_mas(:, :, fr + [-1, 0, 1]);
        I_frames = I_mas(:, :, fr + [-1, 0, 1]);
        Q_frames = Q_mas(:, :, fr + [-1, 0, 1]);
    
        % Identify displacement time t corresponding to current frame fr
        t = min(1 + (fr - 1) * fr2t, length(rea_t) - 1);
        t_0 = floor(t);
        
        % Interpolate consecutive frames in time dimention
        y_mid = find(rea_y >= 0, 1, 'first');

        u_zi = squeeze(u_z(t_0, :, y_mid, :));
        u_zf = squeeze(u_z(t_0 + 1, :, y_mid, :));
    
        t_frac = t - t_0;
        u_rea = (1 - t_frac) * u_zi + t_frac * u_zf;
    
        % Save beamformed frame and true displacement
        save(sprintf('../resources/TrainData/snr%d/ct%4.2f.mat', ...
        img_param.snr, c_t), 'fr', 'RF_frames', 'I_frames', 'Q_frames', ...
        'u_rea', 'rea_x', 'rea_z')

        % Show displacement and center bmode frame
        if graf
        figure(1)
        subplot(1, 2, 1)
        imagesc(est_x, est_z, u_frame)
        colorbar
        subplot(1, 2, 2)
        imagesc(img_x, img_z, I_frames(:, :, 2).^2 + Q_frames(:, :, 2).^2)
        colorbar
        pause()
        end
    end
end