function [TX, TW, TPC] = BMS_Transmit(P, Parameters, Trans)
% Specify TX & TW structure arrays.

% Recover simulation data
lambda = Parameters.lambda;             % Wavelength [m]
n_ang = P.n_ang;                        % Nr. of steering angles
n_push = P.n_push;                      % Number of pushing elements
z_push = P.z_push/ Parameters.lambda;   % Push depth [wvls]
c_push = P.c_push;                      % Nr. of push cycles
n_elem = Trans.numelements;             % Number of elements
elem_space = Trans.spacing;             % Element spacing [wvls]
elem_width = Trans.elementWidth;        % Element width [wvls]
z_foc = Trans.elevationFocusMm * 1e-3;  % Elevation focus depth [m]


% Steered b-mode short pulses

% Calculate steering angles
d_beta = asin(((n_elem-1)*elem_space + elem_width)^-1); % Angle hop [rad]
N_ideal = z_foc / lambda;                               % Ideal # of angles
max_beta = d_beta * N_ideal / 2;                        % Max. angle [rad]

% Specify TX structure array.
TX = repmat(struct( ...
    'waveform', 1, ...                      % selected waveform
    'Origin', [0.0, 0.0, 0.0], ...          % transmit focus origin
    'focus', 0, ...                         % focus depth [wvls]
    'Steer', [0.0, 0.0], ...                % theta, alpha
    'Apod', kaiser(n_elem, 9)', ...         % apodization
    'Delay', zeros(1, n_elem)), ...         % focus delays
    1, n_ang);  

for n = 1:n_ang
    % Calculate current angle
    if n_ang > 1; beta = (2 * n - 1 - n_ang) / (n_ang - 1) * max_beta;
    else; beta = 0; end

    TX(n).Steer = [beta, 0];
    TX(n).Delay = computeTXDelays(TX(n));
end

TW(1) = struct( ...             % specify waveform
    'type', 'parametric', ...   % waveform design method
    'Parameters', ...           %[freq., duty, #cycl, pol]
        [Trans.frequency, .67, 1, 1], ... 
    'sysExtendBL', 0);          % longer than 25 cycles

TPC(1) = struct( ...            % specify TX power controller profile
    'name', 'Imaging', ...
    'hv', 1.6, ...              % initial bipolar voltage
    'maxHighVoltage', 50, ...   % fixed max. voltage limit
    'highVoltageLimit', 40, ... % variable max. voltage limit
    'xmitDuration', 10, ...     % max. transmit duration [usec]
    'inUse', 1 ...
    );

% Push long pulse
push_apod = [zeros(1, (Trans.numelements - n_push) / 2), ...
            ones(1, n_push), ...
            zeros(1, (Trans.numelements - n_push) / 2)];

TX(n_ang + 1) = struct( ...
    'waveform', 2, ...                      % selected waveform
    'Origin', [0.0, 0.0, 0.0], ...          % transmit focus origin
    'focus', z_push, ...                    % focus depth [wvls]
    'Steer', [0.0, 0.0], ...                % theta, alpha
    'Apod',  push_apod, ...                 % apodization
    'Delay', zeros(1, Trans.numelements));  % focus delays

TX(n_ang + 1).Delay = computeTXDelays(TX(2)); % compute delays

TW(2) = struct( ...             % specify waveform
    'type', 'parametric', ...   % waveform design method
    'Parameters', ...           %[freq., duty, #cycl, pol]
        [Trans.frequency, .5, c_push * 2, 1], ... % 960/2 cycles @ 7.6MHz = 64 us
    'sysExtendBL', 1);          % longer than 25 cycles

TPC(5) = struct( ...            % specify TX power controller profile
    'name', 'Push', ...
    'hv', P.hv, ...             % initial bipolar voltage
    'maxHighVoltage', 90, ...   % fixed max. voltage limit
    'highVoltageLimit', 90, ... % variable max. voltage limit
    'xmitDuration', 200,...     % max. transmit duration [usec]
    'inUse', 0 ...
    );

end
