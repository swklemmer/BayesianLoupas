% This script provides an environment to simulate B-mode images using angle
% compunding and non-linear beamforming

clear all
addpath('./lib/Auxiliar');

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

% Simulate for different material properties
ct_list = 0.25:0.25:3;  % shear wave speed [m/s]
N_exp = 20;              % nr. of experiment per sws

for c_t = ct_list

    % Define Trans structure array
    [P, Trans] = BMS_Trans(P);
    
    % Define system parameters
    Parameters = BMS_Parameters(P, Trans);

    % Define Displacement structure array (simulation only)
    Disp = BMS_Disp(sprintf('./FemData/u_%.0f.h5', 3e3 * c_t^2), Parameters);
    
    for n_i = 1:N_exp
    
        clearvars -except c_t ct_list n_i N_exp P Trans Parameters Disp
        clear functions

        % Define PData structure array
        PData = BMS_PData(P, Trans);
    
        % Redefine Media object at random
        Media = BMS_Media(Parameters, PData);
        
        % Specify Receive Buffers
        [RcvBuffer, InterBuffer] = BMS_Buffers(P, Parameters);
        
        % Gather Parameters, Buffers & DisplaWindows in Resource structure
        Resource = struct( ...
            'Parameters',    Parameters, ...
            'RcvBuffer',     RcvBuffer, ...
            'InterBuffer',   InterBuffer, ...
            'HIFU', struct('voltageTrackP5', 0));
    
        % Specify structure arrays for transmision
        [TX, TW, TPC] = BMS_Transmit(P, Parameters, Trans);
        
        % Specify structure arrays for reception
        [Receive, TGC] = BMS_Receive(P, Trans);
        
        % Specify structure arrays for reconstruction
        [Recon, ReconInfo] = BMS_Recon(P);
        
        % Specify Process structure array
        Process = BMS_Process();
        
        % Specify structure arrays for the Event Sequence
        [SeqControl, Event] = BMS_Event(P);
    
        % Save all the structures to a .mat file and run
        filename = 'BMS';
        save(['MatFiles/', filename]);
        tic();
        VSX;
        elap_t = toc();

        % Save results
        save(sprintf('./BModeData/e%d/ct%4.2f_%d.mat', ...
            P.bmode_dly, c_t, n_i), ...
            'RcvData', 'IData', 'QData', ...
            'img_x', 'img_z', 'BFData', 'elap_t')

        fprintf("c_t: %4.2f, n: %d, t = %3.1f\n", c_t, n_i, elap_t)
    end
end

%% Show sonograms

% Simulate for different material properties
ct_list = 0.25:0.25:3;  % shear wave speed [m/s]
N_exp = 1;              % nr. of experiment per sws
hObject.Value = 1;

for c_t = flip(ct_list)    
    for n_i = 1:N_exp
        % Load results
        load(sprintf('./BModeData/e%d/ct%4.2f_%d.mat', ...
            50, c_t, n_i), ...
            'RcvData', 'IData', 'QData', ...
            'img_x', 'img_z', 'BFData', 'elap_t')

        BMS_show_beamf(hObject, 0, c_t);
    end
end

%% Displacement estimation

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

TDE_estimate_u_lou(IData{1}, QData{1});

% Show displacement
hObject.Value = 1;
BMS_show_movie(hObject, 0);