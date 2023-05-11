function u_time = create_gauss_u(img_param, pos_width, t_line)
%CREATE_GAUSS_U Create a gaussian displacement profile for a given ROI.
% INPUT:
% pos_width (4 x 1): axial [wvls] and lateral [lines] center
% coordinates and standard deviation.
% t_line (t_max x 1): Amplitude of pulse over time

% Generate spatial grid
z_dim = 0:(img_param.f_c * img_param.t_s):img_param.z_max; % [wavelengths]
x_dim = 1:img_param.x_max; % [lines]
[z_grid, x_grid] = meshgrid(z_dim, x_dim);

% Retrieve profile dimensions
z_pos = pos_width(1); % [wvls]
z_wid = pos_width(2); % [wvls]
x_pos = pos_width(3); % [lines]
x_wid = pos_width(4);  % [lines]

% Create normalized profile and expand along time dimention
u_gauss = exp(- (z_grid - z_pos).^2 / z_wid^2 ...
              - (x_grid - x_pos).^2 / x_wid^2); %[wvls]
u_time = u_gauss .* reshape(t_line, 1, 1, []);

end
