function [t_dim, x_dim, y_dim, z_dim, u_x, u_y, u_z] = load_u(fem_file)
%LOAD_U Loads displacements from h5 file.

% Read node positions from h5 file
node_pos = h5read(fem_file, "/Mesh/mesh/geometry");

% Read adquisition times from h5 file
fem_info = h5info(fem_file, "/Function/Displacement");

adq_times = struct2cell(fem_info.Datasets);
adq_times = adq_times(1, :);

% Convert times to floats and sort them
[t_dim, sort_ind] = time2float(adq_times);

% Read node z-displacements from h5 file
u = zeros(length(adq_times), size(node_pos, 2), 3);

for i = 1:length(adq_times)
    vec_u = h5read(fem_file, strcat("/Function/Displacement/", ...
                    adq_times(i)));
    u(i, :, :) = vec_u';
end

% Convert listed data into matrix form
[x_dim, y_dim, z_dim, u_x, u_y, u_z] = scatter3matrix(node_pos, ...
                                                      u(sort_ind, :, :));

end

