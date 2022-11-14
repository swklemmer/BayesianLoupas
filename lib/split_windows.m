function [rf_data, img_param] = split_windows(img_param, rf_lines)
%SPLIT_WINDOWS

% Retrieve imaging parameters
f_c = img_param.f_c;
t_s = img_param.t_s;
z_max = img_param.z_max;
K_len = img_param.K_len;
K_hop = img_param.K_hop;
N = img_param.N;

% Define window parameters
M = ceil(K_len / (f_c * t_s));       % window length [smpls]
M_hop = max(floor(K_hop / (f_c * t_s)), 1); % hop size [smpls]
K = floor((z_max - K_len) / K_hop);  % number of windows

% Save parameters to base
img_param.M = M;
img_param.M_hop = M_hop;
img_param.K = K;

% Create windows
rf_data = zeros(K, M, N);
for k = 1:K
    M_win = (1:M) + (k - 1) * M_hop;
    rf_data(k, :, :) = rf_lines(M_win, :);
end

end
