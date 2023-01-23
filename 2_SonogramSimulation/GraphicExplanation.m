addpath('../lib/SonoSim/')
close all

%% Imaging parameters
img_param = struct(...
    'f_c',      7e6, ...    % [Hz] 
    'bw',       0.5, ...    % [%]
    't_s',      0, ...      % [s]
    'z_max',    20, ...     % axial depth [lmbds]
    'M',        0, ...      % kernel rf length [smpls]
    'M_hop',    0, ...      % kernel rf hop [smpls]
    'M_max',    0, ...      % axial depth [smpls]
    'N',        1, ...      % ensemble length [frames]
    'SNR',      60, ...     % [dB]
    'SNR_rho',  0);         % [peak to peak]

img_param.t_s = 1 / (4 * img_param.f_c);
img_param.SNR_rho = 10^(img_param.SNR/20);

u_true = 0.1; % [lmbds]

% Create gaussian RF Pulse 
[rf_pulse, t_cell] = create_pulse(img_param);

% Create normalized RF Line
rng(420)
rf_lines = create_rf_line(img_param, rf_pulse, t_cell, u_true, 1);
rf_lines = rf_lines ./ (2 * max(abs(rf_lines), [], 'all'));

%% Plot Dummy B-mode image with known displacement

scat_pos = [1.95, 0.1, 3.6, 1.9];

% Create random RF lines
u_true_list = [-0.2, 0.4, -0.5, 0.1];
nr_lines = 4;

for i = 1:nr_lines
    rng(4200 * i)
    rf_lines = create_rf_line(img_param, rf_pulse, t_cell, u_true_list(i));
    evalin('base', sprintf('rf_lines%d = rf_lines;', i))
end

% Create axial and lateral dimensions
img_z = (0:size(rf_lines, 1)-1) / 4;
img_x = (0:nr_lines-1) / 2;
lat_sp = 0.1;

% Plot lines and inter-frame displacement
fig = figure(1);
fig.Position = [100, 100, 400, 400];

subplot(2, img_param.N, 1)
imagesc(img_x, img_z, ...
    [lat_sp * rf_lines2(:, 1) + rf_lines1(:, 1), ...
    lat_sp * (rf_lines1(:, 1) + rf_lines3(:, 1)) + rf_lines2(:, 1), ...
    lat_sp * (rf_lines2(:, 1) + rf_lines4(:, 1)) + rf_lines3(:, 1), ...
    lat_sp * rf_lines3(:, 1) + rf_lines4(:, 1)], ...
    [-3, 3])
colormap(gca, 'bone')
hold on
scatter(img_x, scat_pos, 100, 'r*')
hold off

title('Frame T', 'Interpreter', 'latex')
ylabel('Axial distance [$\lambda$]', 'Interpreter', 'latex')

subplot(2, img_param.N, 2)
i = 2;
imagesc(img_x, img_z, ...
    [lat_sp * rf_lines2(:, i) + rf_lines1(:, i), ...
    lat_sp * (rf_lines1(:, i) + rf_lines3(:, i)) + rf_lines2(:, i), ...
    lat_sp * (rf_lines2(:, i) + rf_lines4(:, i)) + rf_lines3(:, i), ...
    lat_sp * rf_lines3(:, i) + rf_lines4(:, i)], ...
    [-3, 3])
colormap(gca, 'bone')
hold on
scatter(img_x, scat_pos + u_true_list, 100, 'r*')
hold off

title('Frame T+1', 'Interpreter', 'latex')

subplot(2, img_param.N, 3.5)
h = heatmap(img_x, img_z(1:4:end), repmat(u_true_list(1:nr_lines), 5, 1));
colormap(gca, 'default')

h.XLabel = 'Lateral distance [$\lambda$]';
h.YLabel = 'Axial distance [$\lambda$]';
h.Title = 'Inter-frame displacement';
h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
h.NodeChildren(3).Title.Interpreter = 'latex';
annotation('textarrow', [0.8, 0.3], [0.28, 0.5], 'string', ...
    'Axial displacement [$\lambda$]', 'HeadStyle', 'none', 'LineStyle', ...
    'none', 'HorizontalAlignment','center','TextRotation', 270, ...
    'Interpreter', 'latex');

%% Plot lines versus time

fig = figure(1);
fig.Position = [600, 200, 250, 200];
for i = 1:size(rf_lines, 2)
    plot(rf_lines(:, i) + i * 1, 'b')
    hold on
end

hold off
grid on
xlim([1, 40])
xticks([10, 20, 30, 40])
xlabel('Sample t (1\dots M)', 'FontSize', 12, 'Interpreter','latex')

ylim([0.25, img_param.N + 0.75])
yticks(1:img_param.N)
ylabel('Frame T (1\dots N)', 'FontSize', 12, 'Interpreter','latex')

set(gca, 'Ydir', 'reverse')

%% Plot DTFT

f_dim = -0.5:0.01:0.5;
F_dim = -0.5:0.01:0.5;

rf_dtft = zeros(length(f_dim), length(F_dim));
dtft_matrix = zeros(size(rf_lines));

for f = 1:length(f_dim)
    for F = 1:length(F_dim)

        % Create DTFT matrix
        for t = 1:size(rf_lines, 1)
            for T = 1:size(rf_lines, 2)
                dtft_matrix(t, T) = ...
                exp(-1i * 2 * pi * ((t-1) * f_dim(f) + (T-1) * F_dim(F)));
            end
        end
        rf_dtft(f, F) = sum(rf_lines .* dtft_matrix, 'all');
    end
end

fig = figure(2);
fig.Position = [850, 200, 250, 200];

imagesc(f_dim, F_dim, abs(rf_dtft)')
xlabel('Fast frequency f', 'FontSize', 12, 'Interpreter','latex')
ylabel('Slow frequency F', 'FontSize', 12, 'Interpreter','latex')

%% Plot frequency paths

fig2 = gca;
fig = figure(3);
copyobj(fig2, fig)
fig.Position = [1100, 200, 250, 200];
hold on

for u = -0.3:0.05:0.3
    k_u = -u / (img_param.f_c * img_param.t_s);
    if abs(u - u_true) < 1e-4
        plot(f_dim, k_u * f_dim, 'r--', 'LineWidth', 2)
    else
        plot(f_dim, k_u * f_dim, 'r--', 'LineWidth', 0.25)
    end
end
hold off
