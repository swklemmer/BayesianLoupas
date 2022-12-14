function p_xu = eval_likelihood_NCC_2D(rf_data, u_sol)
%LIKELIHOOD_NCC_3D
% Evaluate the likelihood function of the RF data given a candidate
% solution using the "Normalized Cross-Correlation" method.
% size(rf_data) = [# kernels in z, x; # samples in z, x & t]
% size(u_sol) =   [# kernels in z, x]

% For each kernel, evaluate likelihood of current solution
p_xu = zeros(size(rf_data, [1 2]));

for z = 1:size(rf_data, 1)
    for x = 1:size(rf_data, 2)

        % Combine parallel RF lines
        rf_lines = squeeze(sum(rf_data(z, x, :, :, :, :), 5));
        
        % Transform displacement to sample [4 smpls per wvln]
        smpl_lag = - u_sol(z, x) * 4;

        % Calculate normalized cross-correlation
        cross_corr = ...
            xcorr(rf_lines(:, 1), rf_lines(:, 2), 2, 'normalized');

        % Calculate correlation at current solution
        p_xu(z, x) = ...
            interp1(-2:2, cross_corr, smpl_lag, 'spline', 'extrap');
    end
end

end
