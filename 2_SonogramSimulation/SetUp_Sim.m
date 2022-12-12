addpath('../lib/SonoSim/');
addpath('../lib/SonoSim/BMS_aux/');
addpath('../lib/DispEst/')

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

% Transducer (Verasonics L11-5V) and transmit parameters
load('../resources/L11-5v_Steer.mat', 'PData', 'Trans', 'TX');
lambda = 1540 / (Trans.frequency * 1e6);

% Imaging parameters
img_param = struct(...
    'fr_n',     40, ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation (1 --> 20 kHz)
    'snr',      60, ...     % signal-to-noise-ratio [dB]
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    5, ...      % axial kernel hop    [wvls]
    'x_len',    2, ...      % late. kernel length [wvls]
    'x_hop',    2, ...      % late. kernel hop    [wvls]
    't_len',    3 ...       % temp. kernel length [wvls]
    );

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
        'axi_len', 7, ...       % LOU: Axial window length [wvls]
        'axi_hop', 1, ...       % LOU: Axial window hop [wvls]
        'lat_len', 1, ...       % LOU: Lateral window length [wvls]
        'lat_hop', 1, ...       % LOU: Lateral window hop [wvls]
        'med_sz', 7, ...        % LOU: Median filter size [smpls]
        'ens_len', 3, ...       % LOU: N = Ensemble length
        'cum_sum', 15) ...      % LOU: Moving average size [smpls]
        );    