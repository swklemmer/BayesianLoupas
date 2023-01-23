function rf_lines_c = create_rf_line(...
               img_param, rf_pulse, t_cell, u_true, varargin)
%CREATE_RF_LINE Returns single RF lines for a given axial length [s],
% ensemble length, inter-frame displacement [lambda], RF pulse and SNR [dB]

% Retrieve imaging parameters
f_c = img_param.f_c;
bw = img_param.bw;
t_s = img_param.t_s;
z_max = img_param.z_max;
N = img_param.N;
SNR = img_param.SNR;

% Scat position and intensity
M_max = ceil(z_max / (f_c * t_s));
max_u = ceil(abs((N - 1) * u_true / (f_c * t_s))); % [smpls]
n_scat = ceil(0.4 * (M_max + max_u) * t_s / t_cell); % 15 scat/cell
scat_pos = 1 + (1-sign(u_true)) * max_u / 2 + rand(n_scat, 1) * (M_max-1);
scat_int = max(1 + 0.25 * randn(n_scat, 1), 0);

% Create phase-shifted scatterer distributions
scat_map = zeros(M_max + 2 * max_u, N);

for s = 1:n_scat
    for n = 1:N
        % New scatterer position
        pos = scat_pos(s) + (n-1) * u_true / (f_c * t_s); % [smpls]

        % Integer and sub-sample lags
        i_tau = floor(pos); % [smpls]
        sub_tau = (pos - i_tau) * (f_c * t_s); % [rad]
    
        % Add scatterer to map
        scat_map(i_tau, n) = scat_map(i_tau, n) + ...
            scat_int(s) * exp(-sub_tau * 2i * pi);
    end
end

% RF lines
rf_lines = zeros(M_max + 2 * max_u, N);

for n = 1:N
    % Convolution with RF pulse
    rf_lines(:, n) = real(conv(scat_map(:, n), rf_pulse, 'same'));
end

% Normalize RF data and add noise
rf_lines = rf_lines / std(rf_lines(:, 1), 1);
rf_lines_n = rf_lines + 10^(-SNR/20) * randn(M_max + 2 * max_u, N);

% Design hanning windowed band-pass filter
%n_filter = 3; % n taps
%taps_bp = fir1(n_filter, 2 * f_c * t_s * [1 - bw / 2, 1 + bw / 2]);

% Apply band-pass filter
%rf_lines_f = filter(taps_bp, 1, rf_lines_n, [], 1);

% Crop rf_lines
rf_lines_c = rf_lines_n((1:M_max) + max_u, :);

% Show RF lines
if size(varargin) > 0
    fig = figure(1);
    fig.Position = [100, 400, 300, 200];
    stem(scat_pos / 4, scat_int, 'r*', 'MarkerSize', 10)
    xlabel('Axial position [$\lambda$]', 'Interpreter','latex')
    ylabel('Reflectivity [adim]', 'Interpreter','latex')
    xlim([0, (M_max - 1)/4])
    ylim([0, 1.3])
    grid on
    title('\textbf{Scatterers}', 'Interpreter','latex')

%     fig = figure(2);
%     fig.Position = [300, 600, 300, 200];
%     plot((0:M_max-1) * t_s, rf_lines_c);
%     grid on;
%     xlabel('Time [s]')
%     xlim([0, (M_max - 1) * t_s])
%     ylim(sqrt(n_scat) * [-2, 2])
%     title('RF Lines')

    fig = figure(3);
    fig.Position = [300, 600, 300, 200];
    plot((0:M_max-1) / 4, rf_lines_c);
    grid on;
    xlabel('Axial position [$\lambda$]', 'Interpreter','latex')
    xlim([0, (M_max - 1)/4])
    ylim([-3.1, 3.1])
    title('\textbf{Received signal}', 'Interpreter','latex')
    ylabel('Amplitude [V]', 'Interpreter','latex')

end

end
