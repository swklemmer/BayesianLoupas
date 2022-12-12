SetUp_Sim;

ct_list = 0.25:0.25:3;  % shear wave speed [m/s]

% Pre-allocate displacement statistics
mean_u = zeros(length(ct_list), 1);
max_u = zeros(length(ct_list), 1);

for i = 1:length(ct_list)
    c_t = ct_list(i);

    % Load displacement info from .h5 file
    [~, ~, ~, ~, ~, ~, u_z] = ...
        load_u(sprintf('../resources/FemData/u_%.0f.h5', 3e3 * c_t^2));
    
    % Obtatin mean and max displacement
    mean_u(i) = mean(u_z, 'all') / lambda;
    max_u(i) = max(u_z, [], 'all') / lambda;

    histogram(u_z(:) / lambda, [-0.5, 0.5, 10])
    title('Displacement histogram', 'FontSize', 14)
    ylabel("# Nodes", 'FontSize', 12)
    xlabel("Axial Displacement [\lambda]", 'FontSize', 12)
    grid on
    break
end

%% Show results

fig = figure(1);
fig.Position = [0, 0, 750, 300];
subplot(1, 2, 1)
semilogx(3e3 * ct_list.^2, max_u)
title('Peak displacement', 'FontSize', 14)
xlabel("Young's Modulus [Pa]", 'FontSize', 12)
ylabel("Axial Displacement [\lambda]", 'FontSize', 12)
grid on

subplot(1, 2, 2)
semilogx(3e3 * ct_list.^2, mean_u)
title('Mean displacement', 'FontSize', 14)
xlabel("Young's Modulus [Pa]", 'FontSize', 12)
ylabel("Axial Displacement [\lambda]", 'FontSize', 12)
grid on
