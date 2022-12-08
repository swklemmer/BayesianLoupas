function p_xu = eval_likelihood_ACK_2D(rf_data, u_sol)
%LIKELIHOOD_ACK_3D
% Evaluate the likelihood function of the RF data given a candidate
% solution using the "AutoCorrelation Kernel" method.
% size(rf_data) = [# kernels in z, x; # samples in z, x & t]
% size(u_sol) =   [# kernels in z, x]

% Create autocorrelation dimentions
M = size(rf_data, 3);
N = size(rf_data, 5);
tau_m = (1-M):(M-1);
tau_n = (1-N):(N-1);

% For each kernel, evaluate likelihood of current solution
p_xu = zeros(size(rf_data, [1 2]));

for z = 1:size(rf_data, 1)
    for x = 1:size(rf_data, 2)

        % Take mean over lateral dimention
        rf_lines = squeeze(sum(rf_data(z, x, :, :, :), 4));

        % Create Autocorrelation mask [8 smpls per wvl]
        k_u = - u_sol(z, x) * 8;
        phi = sinc(tau_m' + k_u * tau_n);
    
        % Calculate normalized autocorrelation using Wiener-Khinchin Theorem
        gamma = fftshift(ifft2(abs(fft2(rf_lines, 2*M-1, 2*N-1)).^2));
        gamma = gamma / gamma(M, N)^2;
    
        % Accumulate term-by-term multiplication
        p_xu(z, x) = gamma(:)' * phi(:);
    end
end

end