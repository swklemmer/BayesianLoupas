function [cell_vol, x_fwhm, z_fwhm] = resolution_cell(...
                        img_param, sim_param, trans_param, tw, varargin)
%RESOLUTION_CELL Calculate de -6dB elliptic cylinder volume of a given
% transducer configuration.

% Recover simulation data
c_c = sim_param.speedOfSound;                 % Sound Speed [m/s]
lambda = sim_param.lambda;                    % Wavelength [m]
n_ang = img_param.n_ang;                      % Nr. of steering angles
f_s = img_param.f_s;                          % Sampling Frequency [Hz]
f_c = tw.Parameters(1) * 1e6;                 % Center frequency [Hz]
bw = diff(trans_param.Bandwidth * 1e6) / f_c; % Bandwidth [%]
n_elem = trans_param.numelements;             % Number of elements
elem_space = trans_param.spacing;             % Element spacing [lmbds]
elem_width = trans_param.elementWidth;        % Element width [lmbds]
z_foc = trans_param.elevationFocusMm * 1e-3;  % Elevation focus depth [m]

% Create analysis volume
x_dim = (-2:0.02:2) * 1e-3;
y_dim = (-1:0.02:1) * 1e-3;
z_dim = (0:2^14-1) * c_c / f_s / 2;

% Configurate simulation
simulation_config(img_param, trans_param);

% Transducer configuration
trans_tx = create_transducer(img_param, sim_param, trans_param, tw);
trans_rx = create_transducer(img_param, sim_param, trans_param, tw);

% Use all elements for transmission and receive. Use Kaiser window (b = 3)
% for artifact attenuation 
apod = kaiser(n_elem, 9);
xdc_apodization(trans_rx, 0, apod');
xdc_apodization(trans_tx, 0, apod');

% Create pressure field for each steering angle
xdc_center_focus(trans_rx, [0, 0, 0]);
xdc_center_focus(trans_tx, [0, 0, 0]);

% Calculate steering angles
d_beta = asin(((n_elem-1)*elem_space + elem_width)^-1); % Angle hop [rad]
N_ideal = z_foc / lambda;                               % Ideal # of angles
max_beta = d_beta * N_ideal / 2;                        % Max. angle [rad]

% Pre-allocate field measurements
x_das = zeros(length(x_dim), length(z_dim), max(n_ang, 1));
y_das = zeros(length(y_dim), length(z_dim), max(n_ang, 1));

for n = 1:max(n_ang, 1)

    % Calculate current angle
    if n_ang > 1; beta = (2 * n - 1 - n_ang) / (n_ang - 1) * max_beta;
    else; beta = 0; end

    % Create dummy TX structure array
    TX = struct('Origin', [0, 0, 0], ...
               'focus', 0, ...
               'Steer', [beta, 0], ...
               'Apod', ones(1, n_elem, 1));
    
    % Compute delays (must be zero at center element)
    steer_delay = computeTXDelays(TX) / f_c;
    steer_delay = steer_delay - steer_delay(n_elem/2);

    % Steer plane wave at transmission
    xdc_focus_times(trans_tx, 0, steer_delay);

    % Use DAS beamforming at reception
    xdc_dynamic_focus(trans_rx, 0, 0, 0);
    if n_ang == 0; xdc_focus(trans_rx, 0, [0, 0, z_foc]); end   

    % Calculate pulse-echo field in x direction
    for x = 1:length(x_dim)
        [hhp, ~] = calc_hhp(trans_tx, trans_rx, [x_dim(x), 0, z_foc]);
        x_das(x, 1:length(hhp), n) = hhp;
    end
    
    % Calculate pulse-echo field in y direction
    for y = 1:length(y_dim)
        [hhp, ~] = calc_hhp(trans_tx, trans_rx, [0, y_dim(y), z_foc]);
        y_das(y, 1:length(hhp), n) = hhp;
    end
end

% Apply MAS beamforming between angles
if n_ang > 1
    x_mas = zeros(length(x_dim), length(z_dim));
    y_mas = zeros(length(y_dim), length(z_dim));
    
    for m = 1:(n_ang-1)
        for n = (m+1):n_ang
            x_mas = x_mas + ...
                sign(x_das(:, :, m) .* x_das(:, :, n)) .* ...
                sqrt(abs(x_das(:, :, m) .* x_das(:, :, n)));
    
            y_mas = y_mas + ...
                sign(y_das(:, :, m) .* y_das(:, :, n)) .* ...
                sqrt(abs(y_das(:, :, m) .* y_das(:, :, n)));
        end
    end

    % Apply band-pass filter
    bp_fir = fir1(201, 4 * f_c / f_s * (1 + [-bw, bw] / 2), 'bandpass');
    x_fil = filtfilt(bp_fir, 1, x_mas')';
    y_fil = filtfilt(bp_fir, 1, y_mas')';

    % Demodulate RF signals at 2 f_c
    x_demod = abs(demodulate_rf(2 * f_c, f_s, 201, x_fil'))';
    y_demod = abs(demodulate_rf(2 * f_c, f_s, 201, y_fil'))';

else
    % For a single angle
    x_mas = squeeze(x_das(:, :, 1));
    y_mas = squeeze(y_das(:, :, 1));
    x_fil = x_mas;

    % Demodulate RF signals
    x_demod = abs(demodulate_rf(f_c, f_s, 201, x_mas'))';
    y_demod = abs(demodulate_rf(f_c, f_s, 201, y_mas'))';
end

% Normalize IQ signals
x_demod = x_demod ./ max(x_demod, [], 'all');
y_demod = y_demod ./ max(y_demod, [], 'all');

% Find center lines
z_c = x_demod((length(x_dim)+1)/2, :);
[~, c_line] = find(z_c == max(z_c, [], 'all'), 1, 'first');
x_c = x_demod(:, c_line);
y_c = y_demod(:, c_line);

% Calculate -6dB widths in lateral and elevation direction
x_0 = find(x_c > 10^(-6/20), 1, "first");
x_f = find(x_c > 10^(-6/20), 1, "last");
x_fwhm = x_dim(x_f) - x_dim(x_0);

y_0 = find(y_c > 10^(-6/20), 1, "first");
y_f = find(y_c > 10^(-6/20), 1, "last");
y_fwhm = y_dim(y_f) - y_dim(y_0);

% Calculate -6dB width of transmited pulse
pulse_env = abs(demodulate_rf(f_c, 250e6, 61, tw.Wvfm2Wy));
pulse_env = (pulse_env / max(pulse_env, [], 1));
z_0 = find(pulse_env > 10^(-6/20), 1, "first");
z_f = find(pulse_env > 10^(-6/20), 1, "last");
z_fwhm = z_dim(z_f - z_0);

% Calculate cell volume
cell_vol = pi / 4 * x_fwhm * y_fwhm * z_fwhm;

% Free memomry
xdc_free(trans_tx);
xdc_free(trans_rx);

% Show results
if ~isempty(varargin)

fig = figure(1);
plot(x_dim, x_c);
fig.Position = [800, 600, 500, 350];
yline(10^(-6/20), '--')
xline(x_dim([x_0, x_f]), '--')
ylim([0, 1])
xlabel("Lateral distance [m]", 'FontSize', 14)
ylabel("Normalized transducer response", 'FontSize', 14)
title("-6dB resolution cell in X direction", 'FontSize', 16)
grid on

fig = figure(2);
plot(y_dim, y_c);
fig.Position = [800, 300, 300, 200];
yline(10^(-6/20), '--')
xline(y_dim(y_0), '--')
xline(y_dim(y_f), '--')
ylim([0, 1.2])
xlabel("Elevation distance [m]")
ylabel("Normalized transducer response")
title("-6dB resolution cell in Y direction")
grid on

fig = figure(3);
plot(z_dim, z_c);
fig.Position = [800, 0, 300, 200];
yline(10^(-6/20), '--')
xline(z_dim(z_0), '--')
xline(z_dim(z_f), '--')
xlim([0, 1.5] * 1e-3)
ylim([0, 1.2])
xlabel("Axial distance [m]")
ylabel("Normalized transducer response")
title("-6dB resolution cell in Z direction")
grid on

fig = figure(4);
plot(z_dim(1:length(pulse_env)), pulse_env);
fig.Position = [0, 600, 300, 200];
yline(10^(-6/20), '--')
xline(z_dim(z_0), '--')
xline(z_dim(z_f), '--')
xlabel("Axial distance [m]")
ylabel("Normalized excitation pulse")
title("-6dB pulse length")
grid on

fig = figure(5);
fig.Position = [0, 300, 300, 200];
imagesc(z_dim, x_dim, x_mas)
xlim([0, 1e-3])
xlabel("Axial distance [m]")
ylabel("Lateral distance [m]")
title("Beamformed RF data")

fig = figure(6);
fig.Position = [0, 0, 300, 200];
imagesc(z_dim, x_dim, x_demod)
xlim([0, 1e-3])
xlabel("Axial distance [m]")
ylabel("Lateral distance [m]")
title("Abs IQ data")

fig = figure(7);
fig.Position = [400, 0, 400, 200];
sgtitle('Frequency Spectrum at center line')
f_dim = linspace(0, 2, length(z_dim));

subplot(1, 3, 1)
title('Beamformed RF data')
plot(f_dim, abs(fft(x_mas((length(x_dim)+1)/2, :))))
xline(4 * f_c / f_s * (1 + [-bw, bw] / 2))
xlim([0, 0.2])
ylim([0, 2e-23])

subplot(1, 3, 2)
title('Filtered RF data')
plot(f_dim, abs(fft(x_fil((length(x_dim)+1)/2, :))))
xline(4 * f_c / f_s * (1 + [-bw, bw] / 2))
xlim([0, 0.2])
ylim([0, 2e-23])

subplot(1, 3, 3)
title('Demodulated RF data')
plot(f_dim, abs(fft(x_demod((length(x_dim)+1)/2, :))))
xline(bw / 2 * f_c / f_s)
xlim([0, 0.2])
ylim([0, 2e-23])
end
end

