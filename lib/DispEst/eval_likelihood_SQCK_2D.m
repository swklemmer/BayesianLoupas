function p_xu = eval_likelihood_SQCK_2D(img_p, rf_data, u_sol)
%LIKELIHOOD_SQCK_3D
% Evaluate the likelihood function of the RF data given a candidate
% solution using the "Search Quick AutoCorrelation Kernel" method.
% size(rf_data) = [# kernels in z, x; # samples in z, x & t]
% size(u_sol) =   [# kernels in z, x]

persistent max_fc f_c_f_s

% Retrieve maximum fast frequency column
param_flag = evalin('base', 'param_flag');
if param_flag
    evalin('base', 'param_flag = 0;')
    max_fc = evalin('base', 'max_fc');
    f_c_f_s = img_p.f_c * img_p.t_s;
end

% Create autocorrelation dimentions
M = size(rf_data, 4);
N = size(rf_data, 6);
tau_m = 0:(M-1);
tau_n = 0:(N-1);

% For each kernel, evaluate likelihood of current solution
p_xu = zeros(size(rf_data, [1 2]));

for z = 1:size(rf_data, 1)
    for x = 1:size(rf_data, 2)
        % Combine parallel RF lines
        rf_lines = squeeze(sum(rf_data(z, x, :, :, :, :), 5));

        % Create Autocorrelation mask 
        phi = exp(-2j * pi * max_fc(z, x) * (tau_m' - u_sol(z, x) / f_c_f_s * tau_n));
    
        % Evaluate DTFT at desired point
        p_xu(z, x) = abs(rf_lines(:)' * phi(:)) / (M * N);
    end
end

end
