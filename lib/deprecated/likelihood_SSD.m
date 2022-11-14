function [p_xu, elapsed_t] = likelihood_SSD(...
                f_c, t_s, M, N, SNR_rho, rf_lines, u_dim, alpha, varargin)
%LIKELIHOOD_SQD
% Returns likelihood function using SSD (Sum of Squared Differences).
% Uses spline interpolation to shift signals with sub-sample resolution.
% When N > 2, the results of various frame pairs are averaged.


% Run only when parameters are changed
if evalin('base', 'param_flag')

    % Lower parameter change flag
    assignin('base', 'param_flag', 0);
end

% Calculate signal power
snr_rho = 10^(0.1/20);
lambda = 7e0;
sigma = mean(rms(rf_lines).^2) / (1 + lambda);

tic()

% Calculate likelihood
p_xu = zeros(length(u_dim), 1);

for n = 1:(N-1)
    for u = 1:length(u_dim)
        % Translate second frame in u [lmbds] using spline interpolation
        rf_trans = fraccircshift(rf_lines(:, n+1), - u_dim(u) / (f_c * t_s));

        % Accumulate square differences
        p_xu(u) = p_xu(u) + ...
            exp(- mean((rf_lines(:, n) - rf_trans).^2) / (2 * sigma));
    end
end

% Normalize likelihood
p_xu = p_xu / (N - 1) / (2 * pi * sigma)^(M/2);

% Return elapsed time
elapsed_t = toc();

% Figures
if ~isempty(varargin{1})

    u_true = cell2mat(varargin{1});

    fig = figure(9);
    fig.Position = [600, 350, 300, 200];

    plot(u_dim, p_xu)
    xline(u_true, '--')
    xlabel('Disp. [lambda]')
    ylim([0, 4])
    title('SSD: Likelihood')
    grid on
end

end
