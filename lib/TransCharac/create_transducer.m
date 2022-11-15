function trans_obj = create_transducer(img_param, sim_param, trans, tw)
%STANDARD_TRANSDUCER Creates a Field II Transducer object based on a list
% of transducer characteristics.

n_elem = trans.numelements;                     % Number of elements
width = trans.elementWidth * sim_param.lambda;  % Element width [m]
kerf = trans.spacingMm * 1e-3 - width;          % Element spacing [m]
height = trans.elevationApertureMm * 1e-3;      % Element height [m]
elv_focus = trans.elevationFocusMm * 1e-3;      % Elevation Focus [m]

% Define transducer object
trans_obj = xdc_focused_array(n_elem, width, height, kerf, elv_focus, ...
                          1, 10, [0 0 0]);

% Sampling and central frecuency
f0 = trans.frequency * 1e6; % [Hz]
fs = img_param.f_s;         % [Hz]

% Retrieve transducer's impulse response
IR = trans.IR1wy;
t_IR = (0:(length(IR)-1)) / 250e6;

% Interpolate IR at sampling frequency
t_ip = 0:1/fs:4.5/f0;
IR_ip = interp1(t_IR, IR, t_ip, 'linear', 0);
xdc_impulse(trans_obj, IR_ip)

% Retrieve transducer's input voltage signal
TW = tw.TriLvlWvfm_Sim;
t_TW = (0:(length(TW)-1)) / 250e6;

% Interpolate voltage signal at sampling frequency
t_ip = 0:1/fs:2/f0;
TW_ip = interp1(t_TW, TW, t_ip, 'nearest', 0);
xdc_excitation(trans_obj, TW_ip);

% Activate baffle
xdc_baffle(trans_obj, 1);
end

