function rf_batch_n = create_rf_batch(img_param, u_true, varargin)
%CREATE_RF_LINE Returns a batch of RF lines for a given axial length [s],
% ensemble length, inter-frame displacement [lambda], RF pulse and SNR [dB]

% Retrieve imaging parameters
f_c = img_param.f_c;
t_s = img_param.t_s;
z_max = ceil(img_param.z_max / (f_c * t_s)); % [smpls]
x_max = img_param.x_max;
t_max = img_param.t_max;
snr = img_param.snr;

% Create pulse
[rf_pulse, t_cell] = create_pulse(img_param);

% Calculate output array dimentions
max_u = ceil(max(u_true, [], 'all') / (f_c * t_s)); % [smpls]
max_gap = (1 + sign(max_u)) / 2 * abs(max_u); % [smpls]
min_u = floor(min(u_true, [], 'all') / (f_c * t_s)); % [smpls]
min_gap = (1 - sign(min_u)) / 2 * abs(min_u); % [smpls]

% Scatterer number
n_scat = ceil(15 * z_max * t_s / t_cell); % 15 scat/cell/line

% Create phase-shifted scatterer distributions
scat_map = zeros(min_gap + z_max + max_gap, x_max, t_max);

for x = 1:x_max
    % Scatterer z-position
    scat_pos = 1 + min_gap + rand(n_scat, 1) * z_max;
    
    % Scaterrer intensity
    scat_int = max(1 + 0.25 * randn(n_scat, 1), 1e-3);

    for s = 1:n_scat
        for t = 1:t_max
    
            % Interpolate displacement at scat. position
            u_inter = interp2(0:(z_max-1), 1:x_max, u_true(:, :, t), ...
                scat_pos(s) - min_gap, x, 'linear', 0);
    
            % New scatterer position
            pos = scat_pos(s) + u_inter / (f_c * t_s); % [smpls]
    
            % Integer and sub-sample lags
            i_tau = floor(pos); % [smpls]
            sub_tau = (pos - i_tau) * (f_c * t_s); % [rad]
        
            % Add scatterer to map
            scat_map(i_tau, x, t) = scat_map(i_tau, x, t) + ...
                                    scat_int(s) * exp(-sub_tau * 2i * pi);
        end
    end
end

% RF batch
rf_batch = zeros(min_gap + z_max + max_gap, x_max, t_max);

for t = 1:t_max
    for x = 1:x_max
        % Convolution with RF pulse
        rf_batch(:, x, t) = real(conv(scat_map(:, x, t), rf_pulse, 'same'));
    end
end

% Crop data
rf_batch_c = rf_batch(min_gap + (1:z_max), :, :);

% Normalize RF data and add noise
rf_batch_c = rf_batch_c / std(rf_batch_c(:, 1, 1), 1);
rf_batch_n = rf_batch_c + 10^(-snr/20) * randn(z_max, x_max, t_max);

% Show RF lines
if size(varargin) > 0
    fig = figure(1);
    fig.Position = [100, 400, 300, 200];
    stem(scat_pos / 4, scat_int, 'r*', 'MarkerSize', 10)
    xlabel('Axial position [$\lambda$]', 'Interpreter','latex')
    ylabel('Reflectivity [adim]', 'Interpreter','latex')
    xlim([0, (z_max - 1)/4])
    ylim([0, 1.3])
    grid on
    title('\textbf{Scatterers}', 'Interpreter','latex')

    fig = figure(2);
    fig.Position = [300, 600, 300, 200];
    plot((0:z_max-1) * t_s, rf_batch_n);
    grid on;
    xlabel('Time [s]')
    xlim([0, (z_max - 1) * t_s])
    ylim(sqrt(n_scat) * [-2, 2])
    title('RF Lines')

    fig = figure(3);
    fig.Position = [300, 600, 300, 200];
    plot((0:z_max-1) / 4, rf_batch_n);
    grid on;
    xlabel('Axial position [$\lambda$]', 'Interpreter','latex')
    xlim([0, (z_max - 1)/4])
    ylim([-3.1, 3.1])
    title('\textbf{Received signal}', 'Interpreter','latex')
    ylabel('Amplitude [V]', 'Interpreter','latex')

end

end
