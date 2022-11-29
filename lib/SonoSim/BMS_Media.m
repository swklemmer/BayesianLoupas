function Media = BMS_Media(Parameters, PData)

x_span = PData(1).Size(2) * PData(1).PDelta(1);     % x-span [wvls]
y_span = 1.04e-3 / Parameters.lambda;               % y-span [wvls]
z_span = PData(1).Size(1) * PData(1).PDelta(3);     % z-span [wvls]
x_start = PData(1).Origin(1);                       % x starting pos.  [wvls]
z_start = PData(1).Origin(3);             % z starting depth [wvls]

% Number of scaterrers (density = 10 scat./cell)
res_cell = 0.027 / (Parameters.lambda * 1e3)^3; % wvls^3
num_points = ceil(x_span * z_span * y_span / res_cell * 10);
Media.numPoints = num_points;

% Place scatteres at random
Media.MP = zeros(num_points, 4);
Media.MP(:, 1) = x_start + rand(num_points, 1) * x_span; % x-position
Media.MP(:, 2) = (rand(num_points, 1) - 0.5) * y_span; % y-position
Media.MP(:, 3) = z_start + rand(num_points, 1) * z_span; % z-position
Media.MP(:, 4) = max(0.5 + rand(num_points, 1), 1e-3);   % reflectivity

Media.OP = Media.MP;              % original position
Media.attenuation = -0.3;         % attenuation in dB/MHz/cm
Media.function = 'BMS_move_scat'; % scatterer movement function

end
