function [x_dim, y_dim, z_dim, u_x, u_y, u_z] = scatter3matrix(node_pos, ...
                                                          u_scatter)
%SCATTER3MATRIX Transforms displacement in scatter vector form to 4D matrix
%representation.

x_dim = unique(node_pos(1, :));
y_dim = unique(node_pos(2, :));
z_dim = unique(node_pos(3, :));

dx = x_dim(2) - x_dim(1);
dy = y_dim(2) - y_dim(1);
dz = z_dim(2) - z_dim(1);

u_x = zeros(size(u_scatter, 1), length(x_dim), length(y_dim), length(z_dim));
u_y = zeros(size(u_scatter, 1), length(x_dim), length(y_dim), length(z_dim));
u_z = zeros(size(u_scatter, 1), length(x_dim), length(y_dim), length(z_dim));

for t = 1:size(u_scatter, 1)
    for i = 1:length(u_scatter)

        % Find position indices
        x_ind = round(1 + node_pos(1, i) / dx);
        y_ind = round(1 + node_pos(2, i) / dy);
        z_ind = round(1 + node_pos(3, i) / dz);

        % Save displacement
        u_x(t, x_ind, y_ind, z_ind) = u_scatter(t, i, 1);
        u_y(t, x_ind, y_ind, z_ind) = u_scatter(t, i, 2);
        u_z(t, x_ind, y_ind, z_ind) = u_scatter(t, i, 3);

    end
end
end

