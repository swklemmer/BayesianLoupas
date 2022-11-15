function simulation_config(img_param, trans)
%SIMULATION_CONFIG Sets up Field II with following parameters
% Returns: Acoustic Impedance and Absorption

% Recover parameters
f0 = trans.frequency;         % Transducer center frequency [Hz]
fs = img_param.f_s;           % Sampling frequency [Hz]
alpha = img_param.att * 1e-4; % Attenuation Coeficient [dB/(m*Hz)]

% Set Field II's sampling frequency
set_sampling(fs);

% Configurate Attenuation
set_field('att_f0', f0)         % Absorption center frecuency
set_field('freq_att', alpha)    % [dB/(cm * MHz)]
set_field('att', f0 * alpha)    % Frequency dependent attenuation
set_field('use_att', 1);

end
