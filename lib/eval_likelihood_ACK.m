function p_xu = eval_likelihood_ACK(img_param, rf_data, u_sol)
%LIKELIHOOD_ACK 
% size(rf_data) = K, N, M
% size(u_sol) = K, 1

% Retrieve imaging parameters
f_c = img_param.f_c;
t_s = img_param.t_s;
K = img_param.K;
M = img_param.M;
N = img_param.N;
%SNR_rho = img_param.SNR_rho;

% For each kernel, evaluate likelihood of current solution
p_xu = zeros(K, 1);

tau_m = (1-M):(M-1);
tau_n = (1-N):(N-1);

for k = 1:K

    % Create Autocorrelation mask
    k_u = - u_sol(k) / (f_c * t_s);
    phi = sinc(tau_m' + k_u * tau_n);

    % Calculate normalized autocorrelation using Wiener-Khinchin Theorem
    gamma = fftshift(ifft2(abs(fft2(squeeze(rf_data(k, :, :)),...
        2*M-1, 2*N-1)).^2));
    gamma = gamma / gamma(M, N)^2;

    p_xu(k) = gamma(:)' * phi(:);
end

end