function demod_data = demod_batch(rf_data, f_c, f_s)
%DEMOD_BATCH Demodulates RF data from Telemed acquisitions.
% INPUT: rf_data (size: z_len * x_len * t_len)
%        freqUS  (Center frequency)
%        FR      (Sampling frequency)
% OUTPUT: demod_data (size: z_len * x_len * t_len)


% Retrieve signal dimentions
z_len = size(rf_data, 1);
x_len = size(rf_data, 2);
t_len = size(rf_data, 3);

% Mix signal down using local oscilator
z_line = (0:z_len-1) / f_s; % fast time [s]
down_mix = 2 * rf_data .* ...
    repmat(exp(-2i * pi * f_c * z_line'), 1, x_len, t_len);

% Design hanning windowed low-pass filter
taps_lp = fir1(ceil(z_len / 10), 1.7 * f_c / f_s, 'low');

% Apply low-pass filter
demod_data = filtfilt(taps_lp, 1, down_mix);
end

