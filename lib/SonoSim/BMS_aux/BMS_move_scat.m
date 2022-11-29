function BMS_move_scat

% Declare persistent variables to track time
persistent fr fr_max t_dim fr2t x_dim y_dim z_dim

% Run only at initialization
if isempty(fr)
    fr = 1;
    fr_max = evalin('base', 'P.bmode_adq');

    % Get displacement dimentions
    t_dim = evalin('base', 'Disp.t_dim');
    x_dim = evalin('base', 'Disp.x_dim');
    y_dim = evalin('base', 'Disp.y_dim');
    z_dim = evalin('base', 'Disp.z_dim');

    % Calculate frame rate to disp. data time resolution ratio
    fr2t = evalin('base', 'P.bmode_dly') * 1e-6 / diff(t_dim([1, 2]));
end


% Identify data index t corresponding to current frame fr
t = min(1 + fr * fr2t, length(t_dim) - 1);
t_0 = floor(t);

% Retrieve scatterer position and displacement info
Media = evalin('caller', 'Media');
u_xi = squeeze(evalin('base', sprintf('Disp.u_x(%d, :, :, :)', t_0)));
u_xf = squeeze(evalin('base', sprintf('Disp.u_x(%d, :, :, :)', t_0 + 1)));
u_yi = squeeze(evalin('base', sprintf('Disp.u_y(%d, :, :, :)', t_0)));
u_yf = squeeze(evalin('base', sprintf('Disp.u_y(%d, :, :, :)', t_0 + 1)));
u_zi = squeeze(evalin('base', sprintf('Disp.u_z(%d, :, :, :)', t_0)));
u_zf = squeeze(evalin('base', sprintf('Disp.u_z(%d, :, :, :)', t_0 + 1)));

% Interpolate consecutive frames in time dimention
t_frac = t - t_0;
u_x = (1 - t_frac) * u_xi + t_frac * u_xf;
u_y = (1 - t_frac) * u_yi + t_frac * u_yf;
u_z = (1 - t_frac) * u_zi + t_frac * u_zf;

% Interpolate values in space dimention
Media.MP(:, 1) = Media.OP(:, 1) + ...
        interp3(y_dim, x_dim, z_dim, u_x, abs(Media.OP(:, 2)), ...
        abs(Media.OP(:, 1)), Media.OP(:, 3), 'linear', 0);

Media.MP(:, 2) = Media.OP(:, 2) + ...
        interp3(y_dim, x_dim, z_dim, u_y, abs(Media.OP(:, 2)), ...
        abs(Media.OP(:, 1)), Media.OP(:, 3), 'linear', 0);

Media.MP(:, 3) = Media.OP(:, 3) + ...
        interp3(y_dim, x_dim, z_dim, u_z, abs(Media.OP(:, 2)), ...
        abs(Media.OP(:, 1)), Media.OP(:, 3), 'linear', 0);

% Assign new position values to Media
assignin('base', 'Media', Media);

% Check if end was reached and update time
if fr >= fr_max
    fr = 1;
else
    fr = fr + 1;
end

end
