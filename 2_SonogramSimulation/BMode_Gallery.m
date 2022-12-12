SetUp_Sim;

ct_list = 0.5:0.25:3;  % shear wave speed [m/s]
N_exp = 1;              % nr. of experiment per sws

%% NL-Beamforming over Magnitude Data

for c_t = flip(ct_list)    
    for n_i = 1:N_exp
        % Load results
        load(sprintf('../resources/BModeData/MAS/ct%4.2f_%d.mat', ...
            c_t, n_i), 'img_x', 'img_z', 'BFData')

        % Show displacement
        BMS_show_beamf(hObject, 0, c_t);
    end
end

%% NL-Beamforming over RF Data

for c_t = flip(ct_list)    
    for n_i = 1:N_exp
        % Load results
        load(sprintf('../resources/BmodeData/MAS/ct%4.2f_%d.mat', ...
            c_t, n_i), 'img_x', 'img_z', 'RcvData')

        % Beamform IQ data
        [RF_mas, I_mas, Q_mas] = ...
            beamform_rf_lin(img_param, PData, Trans, TX, RcvData{1});

        % Show displacement
        BFData = I_mas.^2 + Q_mas.^2;
        BMS_show_beamf(hObject, 0, c_t);
    end
end

%% NL-Beamforming over IQ Data

for c_t = flip(ct_list)    
    for n_i = 1:N_exp
        % Load results
        load(sprintf('../resources/BModeData/SPW/ct%4.2f_%d.mat', ...
            c_t, n_i), 'img_x', 'img_z', 'IData', 'QData')

        % Show displacement
        BFData = squeeze(IData{1}.^2 + QData{1}.^2);
        BMS_show_beamf(hObject, 0, c_t);
    end
end