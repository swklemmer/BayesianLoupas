function [est_z, est_x] = est_dimensions(est_p, PData, u_0)
%EST_DIMENSIONS Given certain estimation parameters, obtains the z- and
%x-positions of each kernel 

% Define kernel size
z_len = ceil(est_p.z_len / PData.PDelta(3));       % from wvls to smpls
x_len = ceil(est_p.x_len / PData.PDelta(1));       % from wvls to smpls

% Define kernel hop
z_hop = ceil(est_p.z_hop / PData.PDelta(3));       % from wvls to smpls
x_hop = ceil(est_p.x_hop / PData.PDelta(1));       % from wvls to smpls

% Define number of kernels
z_num = size(u_0, 1);
x_num = size(u_0, 2);

% Create estimation dimensions [wvls]
est_z = ((0:z_num-1) * z_hop + z_len/2) * PData.PDelta(3) + PData.Origin(3);
est_x = ((0:x_num-1) * x_hop + x_len/2) * PData.PDelta(1) + PData.Origin(1);

end

