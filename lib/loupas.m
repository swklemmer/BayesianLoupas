function u = loupas(img_param, iq_lines, varargin)
%LOUPAS 1D

% Separate IQ Data
IData = real(iq_lines);
QData = imag(iq_lines);

% Retrieve imaging parameters
K = img_param.K;
M = img_param.M;
N = img_param.N;
M_hop = img_param.M_hop;

% Estimate velocity
u = zeros(K, 1);

% Calculate ensemble windows
win_t_F = (1:N-1);
win_t_f = (1:N);

for z = 1:K
    % Calculate axial windows
    win_z_F = (1:M) + (z - 1) * M_hop;
    win_z_f = (1:M-1) + (z - 1) * M_hop;

    % Calculate frecuency estimates
    F_0 = atan(...
        sum( ...
        (QData(win_z_F, win_t_F) .* ...
        IData(win_z_F, win_t_F + 1) - ...
        IData(win_z_F, win_t_F) .* ...
        QData(win_z_F, win_t_F + 1)), 'all') / ...
        sum( ...
        (IData(win_z_F, win_t_F) .* ...
        IData(win_z_F, win_t_F + 1) + ...
        QData(win_z_F, win_t_F) .* ...
        QData(win_z_F, win_t_F + 1)), 'all')) / (2*pi);

    f_dopp = atan(...
        sum( ...
        (QData(win_z_f, win_t_f) .* ...
        IData(win_z_f + 1, win_t_f) - ...
        IData(win_z_f, win_t_f) .* ...
        QData(win_z_f + 1, win_t_f)), 'all') / ...
        sum( ...
        (IData(win_z_f, win_t_f) .* ...
        IData(win_z_f + 1, win_t_f) + ...
        QData(win_z_f, win_t_f) .* ...
        QData(win_z_f + 1, win_t_f)), 'all')) / (2*pi);

    % Calculate differential displacement [wvls]
    u(z) =  F_0 / (1 + f_dopp);
end

if ~isempty(varargin)
    fig = figure(3);
    fig.Position = [600, 600, 300, 200];
    u_true = varargin{1};

    plot(u);
    yline(u_true, '--')
    grid on;
    ylim([-0.25, 0.25])
    ylabel('Disp. [lambda]')
    xlabel('Window #')
    title('Loupas dislpacement estimation')
    legend({'estim.', 'real'})

    fprintf('RMSE: %5.3f [smpls]\n', rms(u - u_true))
    fprintf('Rel : %5.1f %% \n', 100 * rms(u - u_true) / abs(u_true))
end

end
