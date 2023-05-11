% This script takes the median among a group of NCC estimations to generate
% a ground truth displacement map that takes the effect of the imaging 
% system into consideration.


%% Load estimates

mkdir('../resources/EstData/ground_truth/')
load_dir = 'ncc_ovs1';
ct_list = 0.75:0.25:3; % [m/s]
tic();

for c_t = 1:length(ct_list)
    for it = 1:10
        % Load estimation
        load(sprintf('../resources/EstData/%s/%.2f_%d.mat', ...
            load_dir, ct_list(c_t), it), ...
            'est_p', 'img_p', 'est_x', 'est_z', 'u_hat', 'err', ...
            'a_gain', 'PData')

        % Pre-allocate estimations containter during first iteration
        if (c_t == 1) && (it == 1)
            u_hats = zeros([size(u_hat), 10]);
        end
        
        % Save estimation
        u_hats(:, :, it) = u_hat;

        fprintf('i =%3d  | c_t =%5.2f | t =%4.1f\n', ...
            it, ct_list(c_t), toc());
    end

    % Take mean among all
    u_mean = squeeze(mean(u_hats, 3));

    % Save mean displacement estimate
    save(sprintf('../resources/EstData/ground_truth/%.2f.mat', ...
       ct_list(c_t)), ...
       'u_mean', 'est_p', 'img_p', 'est_x', 'est_z', 'PData')

end
