function p_u = eval_prior_2D(met_param, u_sol)
%EVAL_PRIOR_2D
% Evaluate prior probability of candidate solution u_sol
% size(u_sol) = [# kernels in z, x]

% Retrieve imaging parameters
max_z = size(u_sol, 1);
max_x = size(u_sol, 2);

% Retrieve method parameters
alpha = met_param.alpha;
p = met_param.p;
vcn_z = met_param.vcn_z;
vcn_x = met_param.vcn_x;

% For each kernel, compute prior
p_u = zeros(size(u_sol));

for z = 1:size(u_sol, 1)
    % Obtain axial vecinity
    win_z = max(1, floor(z-vcn_z/2)) : min(max_z, floor(z+vcn_z/2));

    for x = 1:size(u_sol, 2)
        % Obtain lateral vecinity
        win_x = max(1, floor(x-vcn_x/2)) : min(max_x, floor(x+vcn_x/2));

        % Accumulate L-norm of difference
        p_u(z, x) = mean(...
        (sqrt((u_sol(win_z, win_x) - u_sol(z, x)).^2 + 1e-10)).^p, 'all');
%         p_u(z, x) = mean((u_sol(win_z, win_x) - u_sol(z, x)).^2, 'all');
    end
end

% Scale according to lambda parameter
p_u = - p_u / (p * alpha^p);

end
