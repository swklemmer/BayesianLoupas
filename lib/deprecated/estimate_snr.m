function SNR = estimate_snr(rf_lines, img_param, met_param)
%ESTIMATE_SNR
% Estimate local SNR by finding the maximum value of the correlation
% between first and second frame

% Retrieve imaging and method parameters
f_c = img_param.f_c;
t_s = img_param.t_s;
N = img_param.N;
u_dim = met_param.u_dim;

% Transalte lambda dimention to sample dimention
s_dim = - u_dim / (f_c * t_s);

% Create correlation dimension
max_lag = ceil(1 / (2 * f_c * t_s));
corr_dim = -max_lag:max_lag;

% Save all SNR estimations
rho = zeros(N - 1, 1);

for n = 1:(N-1)
    % Calculate cross correlation of two adyacent frames
    cross_corr = xcorr(rf_lines(:, n), rf_lines(:, n+1), max_lag, 'coeff');
    
    % Interpolate at search region potins
    p_xu = interp1(corr_dim, cross_corr, s_dim, 'spline', 0);
    
    % Find maximum value
    rho(n) =  max(p_xu, [], 'all');
end

% Estimate median SNR and clip at 30 dB
med_rho = median(rho);
SNR = min(med_rho / (1 - med_rho), 10^(1.5));

end

