function [est_z, est_x, RF_kern, I_kern, Q_kern] = split_kernels(...
                                    est_param, PData, RFData, IData, QData)
%SPLIT_KERNELS

% Define kernel size
z_len = ceil(est_param.z_len / PData.PDelta(3));       % from wvls to smpls
x_len = ceil(est_param.x_len / PData.PDelta(1));       % from wvls to smpls
t_len = est_param.t_len;                               % already in smpls

% Define kernel hop
z_hop = ceil(est_param.z_hop / PData.PDelta(3));       % from wvls to smpls
x_hop = ceil(est_param.x_hop / PData.PDelta(1));       % from wvls to smpls

% Define number of kernels
z_num = floor((size(RFData, 1) - z_len) / z_hop + 1);
x_num = floor((size(RFData, 2) - x_len) / x_hop + 1);
t_num = floor(size(RFData, 3) - t_len + 1);

% Preallocate kernels
RF_kern = zeros(z_num, x_num, t_num, z_len, x_len, t_len);
I_kern = zeros(z_num, x_num, t_num, z_len, x_len, t_len);
Q_kern = zeros(z_num, x_num, t_num, z_len, x_len, t_len);

for z = 1:z_num
    z_win = (1:z_len) + (z - 1) * z_hop;
    for x = 1:x_num
        x_win = (1:x_len) + (x - 1) * x_hop;
        for t = 1:t_num
            t_win = (1:t_len) + (t - 1);

            % Save kernel into 6D matrix
            RF_kern(z, x, t, :, :, :) = RFData(z_win, x_win, t_win);
            I_kern(z, x, t, :, :, :) = IData(z_win, x_win, t_win);
            Q_kern(z, x, t, :, :, :) = QData(z_win, x_win, t_win);
        end
    end
end

% Create estimation dimensions [wvls]
est_z = ((0:z_num-1) * z_hop + z_len/2) * PData.PDelta(3) + PData.Origin(3);
est_x = ((0:x_num-1) * x_hop + x_len/2) * PData.PDelta(1) + PData.Origin(1);

end
