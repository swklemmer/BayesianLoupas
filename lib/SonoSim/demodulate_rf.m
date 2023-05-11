function iq_lines = demodulate_rf(f_c, f_s, taps, rf_lines)
%DEMOULATE_RF Returns IQ signals for a given RF ensemble.

% Retrieve imaging dimentions
M = size(rf_lines, 1);
N = size(rf_lines, 2);

% Mix signal down using local oscilator
t_line = (0:M-1) / f_s; % fast time [s]
down_mix = 2 * rf_lines .* repmat(exp(-2i * pi * f_c * t_line'), 1, N);

% Design hanning windowed low-pass filter
taps_lp = fir1(taps, 2 * f_c / f_s, 'low');

% Apply low-pass filter
iq_lines = filtfilt(taps_lp, 1, down_mix);

end
