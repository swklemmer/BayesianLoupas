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

% Sampling frequency and excitation parameters
fs = img_param.f_s;             % [Hz]
f0 = tw.Parameters(1) * 1e6;    % [Hz]
duty = tw.Parameters(2) * 100;  % [%]
n_cylc = tw.Parameters(3) / 2;  % [cycles]

% Retrieve transducer's impulse response
IR = trans.IR1wy;
t_IR = (0:(length(IR)-1)) / 250e6;

% Interpolate IR at sampling frequency
t_ip = 0:1/fs:4.5/f0;
IR_ip = interp1(t_IR, IR, t_ip, 'linear', 0);
xdc_impulse(trans_obj, IR_ip)

% Retrieve transducer's input voltage signal
if n_cylc < 10
    pulse_v = tw.TriLvlWvfm_Sim';
else
    % Create pulse train and polarity
    t_pulse = 0:1/fs:n_cylc/f0;
    puls_tr = square(t_pulse * 4 * pi * f0, duty) + 1;
    pol = square(t_pulse * 2 * pi * f0, 49);
    pulse_v = -puls_tr .* pol / 2;
end

xdc_excitation(trans_obj, pulse_v);

% Activate baffle
xdc_baffle(trans_obj, 1);
end

