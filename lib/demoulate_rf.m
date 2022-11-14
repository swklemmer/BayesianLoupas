function iq_lines = demoulate_rf(img_param, rf_lines)
%DEMOULATE_RF Returns IQ signals for a given RF ensemble.

% Retrieve imaging parameters
f_c = img_param.f_c;
t_s = img_param.t_s;
M_max = size(rf_lines, 1);
N = img_param.N;

% Mix signal down using local oscilator
t_line = (0:M_max-1) * t_s; % fast time [s]
down_mix = 2 * rf_lines .* repmat(exp(-2i * pi * f_c * t_line'), 1, N);

% Design hanning windowed low-pass filter
n_lp = 17; % n taps
taps_lp = fir1(n_lp, 2 * f_c * t_s, 'low');

% Apply low-pass filter
iq_lines = filter(taps_lp, 1, down_mix, [], 1);

end