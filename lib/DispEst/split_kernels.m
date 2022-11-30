function kernel_data = split_kernels(img_param, img_x, img_z, BFData)
%SPLIT_KERNELS

% Define kernel size
z_len = ceil(img_param.z_len / img_z(2)); % from wvls to sampls
x_len = ceil(img_param.x_len / img_x(2)); % from wvls to sampls
t_len = img_param.t_len;                  % already in sampls

% Define kernel hop
z_hop = ceil(img_param.z_hop / img_z(2)); % from wvls to sampls
x_hop = ceil(img_param.x_hop / img_x(2)); % from wvls to sampls

% Define number of kernels
z_num = floor((size(BFData, 1) - z_len) / z_hop);
x_num = floor((size(BFData, 2) - x_len) / x_hop);
t_num = floor(size(BFData, 3) - t_len);

% Preallocate kernels
kernel_data = zeros(z_num * x_num * t_num, z_len, x_len);

for z = 1:z_num
    for x = 1:x_num
        for t = 1:t_num
            break
        end
    end
end

end

