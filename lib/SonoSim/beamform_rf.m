function [RF_mas, IData, QData] = beamform_rf(...
                                    param, PData, Trans, TX, RcvData)
%BEAMFORM_RF Processes the RF acquisitions stored in a single RcvBuffer
% using non-linear DAS-FMAS (Moubark et. al, 2019) at different steering
% angles.

% Normalize RF data
RcvData = double(RcvData) ./ double(max(RcvData, [], 'all'));

% Retrieve process parameters
n_ang = 3;                                          % steering angles
fr_in = param.fr_n;                                 % input frame number
smpls_fr = size(RcvData, 1) / fr_in - 128;          % samples per frame
fr_out = floor(fr_in / param.fr_d) - n_ang + 1;     % output frame number
x_elem = Trans.ElementPos(:, 1);                    % element pos. [wls]
sigma_n = rms(RcvData, 'all') * db2mag(-param.snr); % noise st. deviation

% Split Data into acquisitions
RF_acq = zeros(smpls_fr, size(RcvData, 2), fr_in);

for t = 1:fr_in
    % Add white noise
    RF_acq(:, :, t) = RcvData((1:smpls_fr) + (t - 1) * smpls_fr, :) + ...
                      sigma_n * randn(size(RF_acq, [1, 2]));
end

% Create grid
z_dim = PData.Origin(3) + (0:PData.Size(1) - 1) * PData.PDelta(3);
x_dim = PData.Origin(1) + (0:PData.Size(2) - 1) * PData.PDelta(1);

% Delay-and-Sum beamforming
RF_das = zeros(length(z_dim), length(x_dim), fr_in);

for t = 1:fr_in

    % Element transmition delay
    elem_t = TX(1 + mod(t, n_ang)).Delay; % [wvls]

    for z = 1:length(z_dim)
        for x = 1:length(x_dim)

            % Element independent arrival time [1 wln = 4 smpls]
            tau_0 = interp1(x_elem, elem_t, x_dim(x)) * 4;

            for e = 1:size(RcvData, 2)

                % Element dependent echo time (1 wln = 4 smpls)
                tau = sqrt(z_dim(z)^2 + (x_dim(x) - x_elem(e))^2) * 4;

                tau_int = floor(tau_0 + tau);
                tau_frac = tau_0 + tau - tau_int;
                
                % Accumulate delayed signals
                RF_das(z, x, t) = RF_das(z, x, t) + ... 
                    (1 - tau_frac) * RF_acq(tau_int, e, t) + ...
                          tau_frac * RF_acq(tau_int + 1, e, t);
            end
        end
    end
end

% Interpolate data in frequency domain
RF_das = interpft(RF_das, 2 * size(RF_das, 1));

% Design hanning windowed band-pass filter
taps_bp = fir1(125, [0.3 0.7], 'bandpass');

% Multiply-and-Sum beamforming
RF_mas = zeros([size(RF_das, [1 2]), fr_out]);

for t = 1:fr_out
    s_mas = zeros(size(RF_das, [1 2]));
    
    for m = t:(t+n_ang-2)
        for n = (m+1):(t+n_ang-1)

            % Accumulate product of consecutive frames
            s_mas = s_mas + sign(RF_das(:, :, m) .* RF_das(:, :, n)) .* ...
                sqrt(abs(RF_das(:, :, m) .* RF_das(:, :, n)));
        end
    end

    % Save filtered data to buffer
    RF_mas(:, :, t) = filtfilt(taps_bp, 1, s_mas);
end

% Demodulate RF signals to I+Q baseband
RF_demod = RF_mas .* ...
   repmat(exp(-pi / 2i * (0:size(RF_mas, 1)-1)'), [1, size(RF_mas, 2, 3)]);

% Design hanning windowed low-pass filter
taps_lp = fir1(125, 0.5, 'low');

% Apply low-pass filter
IData = filtfilt(taps_lp, 1, real(RF_demod));
QData = filtfilt(taps_lp, 1, imag(RF_demod));

end
