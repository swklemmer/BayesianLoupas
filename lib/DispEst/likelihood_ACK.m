function [p_xu, elapsed_t] = likelihood_ACK(...
                f_c, t_s, M, N, rf_lines, u_dim, alpha, varargin)
%LIKELIHOOD_ACK 
% Returns likelihood function using ACK (AutoCorrelation Kernels).
% It operates by directly evaluating the power spectrum  for f (fast 
% frequency) and F (slow frequency) values in each candidate displacement
% path. It can be proven that integranting the DTFT along a straight
% frecuency path with slope k_u is equivalent to discretely integrate the
% autocorrelation multiplied by sinc(m + k_u * m). To reduce processing
% time, the Wiener-Khinchin Theorem can be used to calculate the
% autocorrelation in frequency domain.

persistent phi

% Run only when parameters are changed
if evalin('base', 'param_flag')

    % Lower parameter change flag
    assignin('base', 'param_flag', 0);

    % Create Phi Kernel
    phi = zeros(length(u_dim), 2*M-1, 2*N-1);
    tau_m = (1-M):(M-1);
    tau_n = (1-N):(N-1);

    for u = 1:length(u_dim)
        k_u = - u_dim(u) / (f_c * t_s);
        phi(u, :, :) = sinc(tau_m' + k_u * tau_n);
    end
end

tic();
% Calculate normalized autocorrelation using Wiener-Khinchin Theorem
gamma = fftshift(ifft2(abs(fft2(rf_lines, 2*M-1, 2*N-1)).^2));
gamma = gamma / gamma(M, N)^2;

% Calculate likelihood
p_xu = zeros(length(u_dim), 1);

for u = 1:length(u_dim)
    phi_u = squeeze(phi(u, :, :));
    p_xu(u) = exp(min(gamma(:)' * phi_u(:) / (alpha * 1e-2), 700));
end

% Normalize likelihood
p_xu = p_xu / sum(p_xu);

% Return elapsed time
elapsed_t = toc();

% Figures
if ~isempty(varargin{1})

    u_true = cell2mat(varargin{1});

    fig = figure(6);
    fig.Position = [600, 150, 300, 200];

    plot(u_dim, p_xu)
    xline(u_true, '--')
    xlabel('Disp. [lambda]')
    ylim([0, max(p_xu, [], 'all')])
    title('ACK: Likelihood')
    grid on
    hold on

    fig = figure(7);
    fig.Position = [900, 150, 300, 200];

    % Find true displacement Phi kernel    
    [~ , u_best] = find(u_dim >= u_true, 1, 'first');
    imagesc(abs(squeeze(phi(u_best, :, :))))
    title('ACK: Kernel Magnitude @ True Disp.')
    colorbar

    fig = figure(8);
    fig.Position = [1200, 150, 300, 200];

    imagesc(gamma)
    title('ACK: Autocorrelation')
    colorbar
end

end
