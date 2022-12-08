% This script provides an environment to simulate B-mode images using angle
% compunding and non-linear beamforming

cd('../../Vantage-4.7.6/')
activate;
rmpath('lib/')
cd('../BayesianLoupas/2_SonogramSimulation/')
addpath('../lib/SonoSim/');
addpath('../lib/SonoSim/BMS_aux/');

%% Generate Sonograms
clear all

% Specify user defined parameters
P = struct(...
    'startDepth',   25, ...     % acq. start depth [wvls]
    'endDepth',     75, ...     % acq. end depth [wvls]
    'latDist',      30, ...     % acq. lateral span [wvls]
    'bmode_dly',    50, ...    % delay between b-mode images [usec] (20kHz)
    'bmode_adq',    2e3 / 50, ... % b-mode acq. number
    'hv',           1.6, ...    % transmition bipolar voltage
    'n_ang',        1, ...      % nr. of steering angles
    'n_push',       32, ...     % nr. of pushing elements
    'z_push',       10e-3, ...  % Push depth [m]
    'c_push',       430, ...    % nr. of push cycles
    'simulate',     1, ...      % enable simulate mode
    'c_t',          1.5 ...     % simulated shear wave speed [m/s]
    );

% Simulate for different material properties
ct_list = 0.5:0.25:3;  % shear wave speed [m/s]
N_exp = 1;              % nr. of experiment per sws

for c_t = ct_list

    % Define Trans structure array
    [P, Trans] = BMS_Trans(P);
    
    % Define system parameters
    Parameters = BMS_Parameters(P, Trans);

    % Define Displacement structure array (simulation only)
    Disp = BMS_Disp(sprintf('../resources/FemData/u_%.0f.h5', ...
        3e3 * c_t^2), Parameters);
    
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
        cd('../../Vantage-4.7.6/')
        save(['MatFiles/', filename]);
        tic(); VSX; elap_t = toc();
        cd('../BayesianLoupas/2_SonogramSimulation/')

        % Save results
        save(sprintf('../resources/BModeData/SPW/ct%4.2f_%d.mat', ...
            c_t, n_i), ...
            'RcvData', 'IData', 'QData', ...
            'img_x', 'img_z', 'BFData', 'elap_t')

        fprintf("c_t: %4.2f, n: %d, t = %3.1f\n", c_t, n_i, elap_t)
    end
end
