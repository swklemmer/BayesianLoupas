addpath('../lib/Field_II/')
addpath('../lib/TransCharac/')
addpath('../lib/')
addpath('../../Vantage-4.7.6/Utilities')
rng(6942069)
field_init(0)


%% Parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('L11-5v_PWI.mat', 'Trans', 'Parameters', 'TW');
Resource.Parameters = Parameters;

% Simulation parameters
img_param = struct(...
    'f_c',      Trans.frequency * 1e6, ...  % Central frequency [Hz]
    'f_s',      250e6, ...                  % Sampling frequency [Hz]
    't_s',      1 / 250e6, ...              % Sampling frequency [Hz]
    'att',      0.3, ...                    % Attenuation [dB/cm/MHz]
    'n_ang',    5, ...                      % Steering angles
    'max_beta', 12, ...                     % Max angle [º]
    'n_push',   64);                        % Push elements

%% Calculate resolution cell volume
param_list = [1:2:13];
val_list = zeros(length(param_list), 1);

for i = 1:length(param_list)
    img_param.n_ang = param_list(i);
    val_list(i) = resolution_cell(img_param, Parameters, Trans, TW);
    fprintf("\n\nRes. Cell. Volume = %.4f mm^3\n\n", val_list(i) * 1e9);
end

%% Show results
figure(1)
plot(param_list, val_list(:, 1) * 1e9, '-*')
ylim([0, 0.05])
xlabel('Max Angle')
ylabel('Lateral Beam Width [mm^3]')
title('Effect of non-linear, multiple angle beamforming')
grid on

%%
figure(2)
plot(param_list, val_list(:, 1) * 1e9, '-*')
ylim([0, 0.05])
xlim([1, 13])
xlabel('# Angles')
ylabel('Lateral Beam Width [mm^3]')
title('Effect of non-linear, multiple angle beamforming')
grid on

%% Terminate program
field_end