function [est_z, est_x] = est_dimensions(est_param, PData, u_size)
%EST_DIMENSIONS Calculates the spatial positions of the samples related to
% a given estimation scenario.

% Define kernel size
z_len = ceil(est_param.z_len / PData.PDelta(3));       % from wvls to smpls
x_len = ceil(est_param.x_len / PData.PDelta(1));       % from wvls to smpls

% Define kernel hop
z_hop = ceil(est_param.z_hop / PData.PDelta(3));       % from wvls to smpls
x_hop = ceil(est_param.x_hop / PData.PDelta(1));       % from wvls to smpls

% Define number of kernels
z_num = u_size(1);
x_num = u_size(2);

% Create estimation dimensions [wvls]
est_z = ((0:z_num-1) * z_hop + z_len/2) * PData.PDelta(3) + PData.Origin(3);
est_x = ((0:x_num-1) * x_hop + x_len/2) * PData.PDelta(1) + PData.Origin(1);

end
