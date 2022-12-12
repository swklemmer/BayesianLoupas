function [RF_comp, IData, QData] = beamform_rf_lin(...
                               param, PData, Trans, TX, RcvData, varargin)
%BEAMFORM_RF Processes the RF acquisitions stored in a single RcvBuffer
% using linear DAS at different steering angles.

% Normalize RF data
RcvData = double(RcvData) ./ double(max(RcvData, [], 'all'));

% Retrieve process parameters
n_ang = param.n_ang;                                % steering angles
fr_in = param.fr_n;                                 % input frame number
smpls_fr = size(RcvData, 1) / fr_in - 128;          % samples per frame
x_elem = Trans.ElementPos(:, 1);                    % element pos. [wls]
sigma_n = rms(RcvData, 'all') * db2mag(-param.snr); % noise st. deviation

if isempty(varargin)
    frames = 1:fr_in;
else
    frames = varargin{1};
end

% Split Data into acquisitions
RF_acq = zeros(smpls_fr, size(RcvData, 2), length(frames));

for t = 1:length(frames)

    % Add white noise
    RF_acq(:, :, t) = RcvData((1:smpls_fr) + (frames(t) - 1) * smpls_fr, :) + ...
                      sigma_n * randn(size(RF_acq, [1, 2]));
end

% Create grid
z_dim = PData.Origin(3) + (0:PData.Size(1) - 1) * PData.PDelta(3);
x_dim = PData.Origin(1) + (0:PData.Size(2) - 1) * PData.PDelta(1);

% Delay-and-Sum beamforming
RF_das = zeros(length(z_dim), length(x_dim), size(RF_acq, 3));

for t = 1:length(frames)

    % Element transmition delay
    elem_t = TX(2 + mod(frames(t), n_ang)).Delay; % [wvls] OJO CON 2 HARDCODEADO

    for z = 1:length(z_dim)
        for x = 1:length(x_dim)

            % Element independent arrival time [1 wln = 4 smpls]
            tau_0 = interp1(x_elem, elem_t, x_dim(x)) * 4;

            for e = 1:size(RcvData, 2)

                % Element dependent echo time (1 wln = 4 smpls)
                tau = sqrt(z_dim(z)^2 + (x_dim(x) - x_elem(e))^2) * 4;

                % Interpolate in time
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

% Apply angle compounding
RF_comp = movmean(RF_das, n_ang, 3, 'Endpoints', 'discard');

% Demodulate RF signals to I+Q baseband
RF_demod = RF_comp .* ...
 repmat(exp(-pi / 2i * (0:size(RF_comp, 1)-1)'), [1, size(RF_comp, 2, 3)]);

% Design hanning windowed low-pass filter
taps_lp = fir1(63, 0.5, 'low');

% Apply low-pass filter
IData = filtfilt(taps_lp, 1, real(RF_demod));
QData = filtfilt(taps_lp, 1, imag(RF_demod));

end
