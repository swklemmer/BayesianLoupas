function [u_tru, mean_u]= interp_real_u(est_x, est_z, max_t, lambda, PData, c_t)
%LOAD_REAL_U

% Load displacement info from .h5 file
[rea_t, rea_x, ~, rea_z, ~, ~, u_z] = ...
    load_u(sprintf('../resources/FemData/u_%.0f.h5', 3e3 * c_t^2));

mean_u = mean(u_z, 'all') / lambda;

% Frame to time transformation
fr2t = 50 * 1e-6 / diff(rea_t([1, 2]));

% Pre-allocate displacement movie data
[x_grid, z_grid] = meshgrid(est_x * lambda, (est_z - PData.Origin(3))* lambda);
u_tru = zeros([size(x_grid), max_t]);

% For each required frame, perform interpolation
for fr = 1:max_t

    % Identify displacement time t corresponding to current frame fr
    t = max(min(1 + (fr - 5) * fr2t, length(rea_t) - 1), 1);
    t_0 = floor(t);
    
    % Interpolate consecutive frames in time dimention
    u_zi = squeeze(u_z(t_0, :, 1, :));
    u_zf = squeeze(u_z(t_0 + 1, :, 1, :));
    
    t_frac = t - t_0;
    u_rea = (1 - t_frac) * u_zi + t_frac * u_zf;
    
    % Interpolate in space dimention
    u_pts = interp2(rea_z, rea_x, u_rea, z_grid(:), x_grid(:), ...
        'linear', 0);
    u_tru(:, :, fr) = reshape(u_pts, size(x_grid)) / lambda;
end
end
