function [p_xu, elapsed_t] = likelihood_NCC(...
                f_c, t_s, rf_lines, u_dim, alpha, varargin)
%LIKELIHOOD_NCC
% Returns likelihood function using NCC (Normalized Cross-Correlation).
% Uses parabolic interpolation to measure likelihood at sumb-sample
% resolution. When N > 2, the results of various frame pairs are averaged.

persistent s_dim max_lag corr_dim

% Run only when parameters are changed
if evalin('base', 'param_flag')

    % Lower parameter change flag
    assignin('base', 'param_flag', 0);

    % Retrieve signal dimensions
    N = size(rf_lines, 2);
    
    % Transalte lambda dimention to sample dimention
    s_dim = - u_dim / (f_c * t_s);
    
    % Create correlation dimension
    max_lag = ceil(1 / (2 * f_c * t_s));
    corr_dim = -max_lag:max_lag;
end

tic()
% Calculate likelihood
p_xu = zeros(length(u_dim), 1);

for n = 1:(N-1)
    % Calculate cross correlation of adyacent frames
    cross_corr = xcorr(rf_lines(:, n), rf_lines(:, n+1),...
        max_lag, 'normalized');

    % Accumulate likelihood functions
    p_xu = p_xu + exp(min(...
        interp1(corr_dim, cross_corr, s_dim, 'spline', 0) / alpha, ...
        700))';
end

% Normalize likelihood
p_xu = p_xu / sum(p_xu) / diff(u_dim(1:2));

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
    ylim([0, 0.1])
    title('NCC: Likelihood')
    grid on
end

end
