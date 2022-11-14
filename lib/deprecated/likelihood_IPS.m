function [p_xu, elapsed_t] = likelihood_IPS(...
                       f_c, t_s, M, N, rf_lines, u_dim, varargin)
%LIKELIHOOD_IPS
% Returns likelihood function using IPS (Interpolated Power Spectrum)
% First, one-sided spectrum is calculated. Then, for each displacement
% canditate, a search path along f (fast freq.) and F (slow freq.) is
% defined. The probability of each displacement is proportional to the sum
% of the interpolated spectrum values along it's corresponding path.
% Chosen interpolation method is 'spline'.

persistent f_dim F_dim f_search F_search

% Run only when parameters are changed
if evalin('base', 'param_flag')

    % Lower parameter change flag
    assignin('base', 'param_flag', 0);
    
    % Define on-sided frequency dimensions
    f_dim = (0:(ceil(M/2) - 1)) / M;
    F_dim = ((0:N-1) - floor(N/2)) / N;
    
    % Define frequency search range
    bw = evalin('base', 'bw');
    f_search = ((-0.75:0.01:0.75) * bw + 1) * f_c * t_s; %[1 / f_s]
    F_search = - u_dim' * f_search / (f_c * t_s);
end

tic();
% One-sided Power Spectrum 
Gamma = abs(fft2(rf_lines)).^2;
Gamma_hat = fftshift(Gamma(1:ceil(M/2), :), 2);

% Calculate likelihood
p_xu = zeros(length(u_dim), 1);

for u = 1:length(u_dim)
    p_xu(u) = sum(interp2(f_dim, F_dim, Gamma_hat',...
        f_search, F_search(u, :), 'linear', 0));
end

% Normalize likelihood
p_xu = p_xu / sum(p_xu);

% Return elapsed time
elapsed_t = toc();

% Figures
if ~isempty(varargin)
    fig = figure(4);
    fig.Position = [0, 150, 300, 200];

    u_true = varargin{1};

    imagesc(f_dim, F_dim, Gamma_hat')
    colorbar
    ax = gca;
    ax.YDir = 'normal';
    ax.XDir = 'normal';
    xlim([0, f_c * t_s * 2])
    ylim([-0.5, 0.5])
    xlabel('f')
    ylabel('F')
    hold on
    
    for u = 1:floor(length(u_dim) / 10):length(u_dim)
        scatter(f_search, F_search(u, :), 'r*')
    end
    [~, u_best] = find(u_dim >= u_true, 1, 'first');
    scatter(f_search, F_search(u_best, :), 'c*')
    hold off
    title('IPS: Frequency search paths')

    fig = figure(5);
    fig.Position = [300, 150, 300, 200];
    plot(u_dim, p_xu)
    xline(u_true, '--')
    xlabel('Disp. [lambda]')
    ylim([0, 0.1])
    title('IPS: Likelihood')
    grid on
end

end
