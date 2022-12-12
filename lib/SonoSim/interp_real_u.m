function [u_true, mean_u]= interp_real_u(est_x, est_z, times, PData, lambda, c_t)
%LOAD_REAL_U Interpolates FEM inter-frame displacement at given spatial and
%temporal positions [in wavelengths].

% Load displacement info from .h5 file
[rea_t, rea_x, ~, rea_z, ~, ~, u_z] = ...
    load_u(sprintf('../resources/FemData/u_%.0f.h5', 3e3 * c_t^2));

mean_u = mean(u_z, 'all') / lambda;

% Frame to time transformation
fr2t = 50 * 1e-6 / diff(rea_t([1, 2]));

% Pre-allocate displacement movie data
[x_grid, z_grid] = meshgrid(est_x * lambda, ...
                           (est_z - PData.Origin(3)) * lambda);
u_true = zeros([size(x_grid), length(times)]);

% For each required frame, perform interpolation
for i = 1:length(times)

    % Identify displacement time t corresponding to current frame fr
    t_1 = max(min(1 + (times(i) - 1) * fr2t, length(rea_t) - 1), 1);
    t_2 = max(min(1 + times(i) * fr2t, length(rea_t) - 1), 1);

    t_1i = floor(t_1);
    t_2i = floor(t_2);
    f1 = t_1 - t_1i;
    f2 = t_2 - t_2i;

    % Interpolate consecutive frames in time dimention
    u_z1i = squeeze(u_z(t_1i, :, 1, :));
    u_z1f = squeeze(u_z(t_1i + 1, :, 1, :));    
    u_z2i = squeeze(u_z(t_2i, :, 1, :));
    u_z2f = squeeze(u_z(t_2i + 1, :, 1, :));
    
    u_intp = (1 - f2) * u_z2i + f2 * u_z2f - (1 - f1) * u_z1i - f1 * u_z1f;
    
    % Interpolate in space dimention
    u_pts = interp2(rea_z, rea_x, u_intp, z_grid(:), x_grid(:), ...
        'linear', 0);
    u_true(:, :, i) = reshape(u_pts, size(x_grid)) / lambda;
end
end
