function [p_xu, elapsed_t] = likelihood_QCK(...
                f_c, t_s, rf_lines, u_dim, alpha, varargin)
%LIKELIHOOD_QCK 
% Returns likelihood function using QCK (Quick AutoCorrelation Kernel).
% Evaluates likelihood at a single point.

% Create DTFT dimentions
tau_m = 0:(size(rf_lines, 1)-1);
tau_n = 0:(size(rf_lines, 2)-1);

% For each kernel, evaluate likelihood of current solution
p_xu = zeros(length(u_dim), 1);

for u = 1:length(u_dim)

    % Create Autocorrelation mask 
    phi_u = exp(-2j * pi * (tau_m' * (f_c * t_s) - u_dim(u) * tau_n));

    % Evaluate DTFT at center frequency
    p_xu(u) = exp(abs(rf_lines(:)' * phi_u(:)) / alpha);
end

tic();

% Normalize likelihood
p_xu = p_xu / sum(p_xu) / diff(u_dim(1:2));

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
    title('QCK: Likelihood')
    grid on
    hold on

end

end
