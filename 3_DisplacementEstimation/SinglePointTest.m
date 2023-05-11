addpath('../../ImageWaves/')
addpath('../lib/DispEst/')
addpath('../lib/SonoSim/')
load('../../resources/AndersMatthias_210601143839.mat')

%% Create synthetic data

% Imaging parameters
img_p = struct(...
    'fr_n',     size(rf, 3), ...     % input frame number
    'fr_d',     1, ...      % frame rate decimation
    'snr',      60, ...     % signal-to-noise-ratio [dB]
    'z_len',    10, ...     % axial kernel length [wvls]
    'z_hop',    5, ...      % axial kernel hop    [wvls]
    'x_len',    1, ...      % late. kernel length [wvls]
    'x_hop',    1, ...      % late. kernel hop    [wvls]
    't_len',    2, ...      % temp. kernel length [wvls]
    'z_max',    size(rf, 1) * freqUS / SR, ...      % axial FOV [wvls]
    'f_c', freqUS, ...
    'bw', 0.8, ...
    't_s', 1/SR, ...
    'N', 2 ...
    );

img_p.M = ceil(img_p.z_max / (img_p.f_c * img_p.t_s));

% Method parameters
met_p = struct(...
    'u_dim', -.5:1e-3:.5, ... % [lmbds]
    'ack_a', 5e-5, ...           % [distribution width]
    'ncc_a', 2e-2, ...           % [distribution width]
    'qck_a', 1e-0);              % [distribution width]

PData = struct(...
    'PDelta', [1, 1, freqUS / (2 * SR)], ...
    'Origin', [0, 0, 0]);

% Establish displacement profile
u_real = (1:size(rf, 2)) / size(rf, 2) - .5;

% Create pulse
[rf_pulse, t_cell] = create_pulse(img_p);

% Create random lines
rng(1234)
rf_synth = zeros(1 + ceil(img_p.z_max / (freqUS / SR)), size(rf, 2), img_p.N);
for i = 1:length(u_real)
    rf_synth(:, i, :) = create_rf_line(img_p, rf_pulse, t_cell, u_real(i));
end

% Demodulate
rf = rf_synth;
iq_lines = reshape(demodulate_rf(freqUS, SR, 100, ...
    squeeze(double(rf(:,:)))), size(rf));

% Split data into kernels
[est_z, est_x, RF_kern, I_kern, Q_kern] = ...
    split_kernels(img_p, PData, rf, real(iq_lines), imag(iq_lines));

%% Observe DTFT and resulting PDF

% Create kernel dimentions
tau_m = 0:(size(RF_kern, 4)-1);
tau_n = 0:(size(RF_kern, 6)-1);

for z = 1:size(RF_kern, 1)
    for x = 1:4:size(RF_kern, 2)
        % Retrieve a single ensembles
        rf_lines = squeeze(sum(RF_kern(z, x, 1, :, :, :), 5));

        % Choose DTFT resolution
        N_dtft = 500;
        f_space = (-N_dtft:2:(N_dtft-1)) / (2 * N_dtft);
        search_f = freqUS / SR + (-.05:1/(size(rf_lines, 1) * 3):.05);

        % Find peak along frequency path
        rf_path = zeros(length(search_f), 1);

        for f = 1:length(search_f)
            % Calculate exponential mask at frequency path
            phi = exp(-2j * pi * (tau_m' * search_f(f) ...
                                + tau_n * search_f(f) * -u_real(x)));
    
            % Evaluate abs(DTFT) at desired point
            rf_path(f) = abs(rf_lines(:)' * phi(:));
        end

        % Find maximum index
        [~, max_f] = max(rf_path);
        max_pdf = find(f_space >= search_f(max_f), 1, 'first');
        
        % Obtain DTFT
        rf_dtft = zeros(N_dtft, N_dtft);
        
        for f = 1:N_dtft
            for v = 1:N_dtft
                % Calculate exponential mask at desired point
                phi = exp(-2j * pi * (tau_m' * ((f - 1) / N_dtft - .5) ...
                                     + tau_n * ((v - 1) / N_dtft - .5)));
        
                % Evaluate abs(DTFT) at desired point
                rf_dtft(f, v) = abs(rf_lines(:)' * phi(:));
            end
        end

        % Show spectrum
        figure(1)
        imagesc(f_space, f_space, rf_dtft')
        line([-(freqUS / SR), (freqUS / SR)], [-1, 1] * -u_real(x), ...
            'LineStyle', '--', 'Color', 'r')
        xline(f_space(max_pdf), 'r--')

        % Show likelihood
        param_flag = 1;
        u_pre = u_real(x);
        [p_xu, elapsed_t] = likelihood(img_p, met_p, 'sqck', rf_lines);

        figure(2)
        plot(met_p.u_dim, p_xu)
        xline(u_real(x), '--')
        xlabel('Disp. [lambda]')
        ylim([0, 40])
        title('SQCK: Likelihood')
        grid on
        hold on

        pause()
    end
end
