addpath('../lib/Field_II/')
addpath('../lib/TransCharac/')
addpath('../lib/')
addpath('../../Vantage-4.7.6/Utilities')
rng(6942069)
field_init(0)

%% Parameters

% Transducer (Verasonics L11-5V) and transmit parameters
load('L11-5v_PWI.mat', 'Trans', 'Parameters', 'TW');

% Simulation parameters
img_param = struct(...
    'f_c',      Trans.frequency * 1e6, ...  % Central frequency [Hz]
    'f_s',      250e6, ...                  % Sampling frequency [Hz]
    't_s',      1 / 250e6, ...              % Sampling frequency [Hz]
    'att',      0.3, ...                    % Attenuation [dB/cm/MHz]
    'kaiser',   3, ...                      % Kaiser window beta
    'c_push',   430, ...                    % Push cycles
    'z_push',   10e-3, ...                  % Push focus depth [m]
    'n_push',   64);                        % Push elements

% Media parameters
c_c = Parameters.speedOfSound;  % Speed of sound [m/s]
rho = 1000;                     % Density [kg/m^3]
imp = rho * c_c;                % Acoustic impedance [Pa s/m^3]
absor = ...                     % Absortion coefficient [Np/m]
    img_param.att * img_param.f_c / 20 * log(10);               

%% Configuration

% Simulation configuration
simulation_config(img_param, Trans);

% Transducer configuration
trans_tx = create_transducer(img_param, Parameters, Trans, TW(2));

%% Acoustic intensity calculation

% Create 0.1 mm grid 
dim_x = (-7:0.1:7) * 1e-3;
dim_z = (0:0.1:25) * 1e-3;

[X, Z] = ndgrid(dim_x, dim_z);
xz_grid = [X(:), zeros(length(X(:)), 1) , Z(:)];

% Focus transducer
xdc_center_focus(trans_tx, [0, 0, 0]);
xdc_focus(trans_tx, 0, [0, 0, img_param.z_push]);

% Choose active elements
n_act = img_param.n_push;
n_elem = Trans.numelements;

apod = [zeros(1, floor((n_elem-n_act)/2)),...
        ones(1, n_act),...
        zeros(1, ceil((n_elem-n_act)/2))];
xdc_apodization(trans_tx, 0, apod);

% Calculate pressure field [Pa = N/m^2]
[pres, t_min] = calc_hp(trans_tx, xz_grid);

% Scale to desired Ispta [W/m^2 = 0.0001 W/cm^2]
Ispta = 1000 * 100^2; % Spatial Peak, Temporal Average Intensity [W/m^2]
avg_inten = reshape(mean(pres.^2 / imp, 1), size(X))';
avg_inten = avg_inten / max(avg_inten, [], 'all') * Ispta;

% Find ARF where intensity is high (> 10% peak)
hi_int_pts = avg_inten > 0.1 * Ispta;
[ind_z, ind_x] = find(hi_int_pts);
f_arf = (2 * absor * avg_inten / c_c) .* hi_int_pts;

%% Exportamos variables para FEniCSx
save('results/ARF.mat', 'f_arf', 'ind_x', 'ind_z', 'dim_x', 'dim_z')

%% Mostramos resultados

load('results/ARF.mat', 'f_arf', 'ind_x', 'ind_z', 'dim_x', 'dim_z')
[X, Z] = ndgrid(dim_x, dim_z);

figure(3)
mesh(X, Z, avg_inten')
xlabel("Distancia lateral [m]")
ylabel("Distancia axial [m]")
title("Intensidad promedio-temporal [W/m^2]")

figure(4)
imagesc(dim_x, dim_z, avg_inten)
xlabel("Distancia lateral [m]")
ylabel("Distancia axial [m]")
title("Intensidad promedio-temporal [W/m^2]")
axis image
colorbar

figure(5)
mesh(X, Z, f_arf')
xlabel("Distancia lateral [m]")
ylabel("Distancia axial [m]")
title("Fuerza de radiación acústica [N/m^3]")

%% Finalizamos
xdc_free(trans_tx)
field_end

