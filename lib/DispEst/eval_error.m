function e_met = eval_error(u_real, u_est, est_z, est_x)
%EVAL_ERROR Calculates the following error metrics between the true
%displacement profile and the estimated one. It takes into consideration
%the underestimation and bulrring effects incorporated by the imaging
%system, as studied by [1]
% [1]  Palmeri et. al (2006). Ultrasonic Tracking of Acoustic Radiation Force-induced Displacements
% in Homogeneous Media

% Pre-allocate error metrics
e_met = zeros(size(u_est, 3), 3);

% Pre grid for interpolation
[z_grid, x_grid] = meshgrid(est_z, est_x);
    
for t = 1:size(u_real, 3)

    % Interpolate real displacement at estimation points
    u_intp = interp2(z_real, x_real, u_real(:, :, t), ...
        z_grid(:), x_grid(:), 'linear', 0);

    u_true = reshape(u_intp, size(u_est, 1:2));

    % Apply blurring according with the system's PSF
    %u_true = u_true;

    % Apply underestimation according to 1/F number
    %u_true = u_true;

    % Save error metrics
    e_met(t, 1) = mean(u_hat - u_true, 'all');
    e_met(t, 2) = var(u_hat - u_true, [], 'all');
    e_met(t, 3) = sqrt(e_met(2) + e_met(1)^2);

end


end

