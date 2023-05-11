function [p_xu, elapsed_t] = likelihood_SQCK(...
                f_c, t_s, rf_lines, u_dim, alpha, varargin)
%LIKELIHOOD_QCK 
% Returns likelihood function using SQCK (Search QCK).
% Evaluates likelihood at a single point, but the fast frequency position
% is selected based on the peak along the initial frequency path.

% Obtain displacement estimate using Loupas
u_pre = evalin('base', 'u_pre');

% Find peak along frequency path
search_f = f_c * t_s + (-.05:1/(size(rf_lines, 1) * 3):.05);
rf_path = zeros(length(search_f), 1);

% Create DTFT dimentions
tau_m = 0:(size(rf_lines, 1) - 1);
tau_n = 0:(size(rf_lines, 2) - 1);

for f = 1:length(search_f)
    % Calculate exponential mask at frequency path
    phi = exp(-2j * pi * (tau_m' * search_f(f) ...
                        + tau_n * search_f(f) * -u_pre));

    % Evaluate abs(DTFT) at desired point
    rf_path(f) = abs(rf_lines(:)' * phi(:));
end

% Find maximum index
[~, max_index] = max(rf_path);
max_fc = search_f(max_index);


% Evaluate likelihood for each candidate displacement
p_xu = zeros(length(u_dim), 1);

for u = 1:length(u_dim)

    % Create Autocorrelation mask 
    phi_u = exp(-2j * pi * max_fc * (tau_m' - u_dim(u) / (f_c * t_s) * tau_n));

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
