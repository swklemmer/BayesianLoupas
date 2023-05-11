function u = loupas_3D(est_p, I_kern, Q_kern)
%LOUPAS_3D Calculate initial displacement estimation for a collection of
% 3D kernels using Loupas' algorithm. The reference frame can either be the
% first frame or a moving reference.

% Preallocate displacement data
u = zeros(size(I_kern, [1, 2, 3]));

% If the reference frame is fixed, save all kernels at first time step
if ~est_p.fin_mref
    I_ref = reshape(I_kern(:, :, 1, :, :, 1), size(I_kern, [1 2 4 5]));
    Q_ref = reshape(Q_kern(:, :, 1, :, :, 1), size(I_kern, [1 2 4 5]));
end

% Calculate axial and ensemble windows
win_z_F = 1:size(I_kern, 4);
win_z_f = 1:(size(I_kern, 4)-1);
win_t_F = 1:(size(I_kern, 6)-1);
win_t_f = 1:size(I_kern, 6);

for z = 1:size(I_kern, 1)
    for x = 1:size(I_kern, 2)
        for t = 1:size(I_kern, 3)

            % Squeeze kernel to 3D
            IData = reshape(I_kern(z, x, t, :, :, :), size(I_kern, 4:6));
            QData = reshape(Q_kern(z, x, t, :, :, :), size(Q_kern, 4:6));

            % If reference is fixed, replace first acquisition in ensemble
            if ~est_p.fin_mref
                IData(:, :, 1) = I_ref(z, x, :, :);
                QData(:, :, 1) = Q_ref(z, x, :, :);
            end

            % Calculate frequency estimates
            F_0 = atan(...
                sum( ...
                (QData(win_z_F, :, win_t_F) .* ...
                IData(win_z_F, :, win_t_F + 1) - ...
                IData(win_z_F, :, win_t_F) .* ...
                QData(win_z_F, :, win_t_F + 1)), 'all') / ...
                sum( ...
                (IData(win_z_F, :, win_t_F) .* ...
                IData(win_z_F, :, win_t_F + 1) + ...
                QData(win_z_F, :, win_t_F) .* ...
                QData(win_z_F, :, win_t_F + 1)), 'all')) / (2*pi);
        
            f_dopp = atan(...
                sum( ...
                (QData(win_z_f, :, win_t_f) .* ...
                IData(win_z_f + 1, :, win_t_f) - ...
                IData(win_z_f, :, win_t_f) .* ...
                QData(win_z_f + 1, :, win_t_f)), 'all') / ...
                sum( ...
                (IData(win_z_f, :, win_t_f) .* ...
                IData(win_z_f + 1, :, win_t_f) + ...
                QData(win_z_f, :, win_t_f) .* ...
                QData(win_z_f + 1, :, win_t_f)), 'all')) / (2*pi);
        
            % Estimate inter-frame displacement [wvls]
            u(z, x, t) = -F_0 / (1 + f_dopp);
        end
    end
end

% If reference is a moving frame, accumulate over time
if est_p.fin_mref
    u = cumsum(u, 3); % [wvls]
end

end
