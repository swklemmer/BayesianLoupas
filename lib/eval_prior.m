function p_u = eval_prior(img_param, met_param, u_sol)
%EVAL_PRIOR 
% Evaluate prior probability of candidate solution u_sol

% Retrieve imaging parameters
K = img_param.K;

% Retrieve method parameters
lambda = met_param.lambda;
p = met_param.p;
B = met_param.B;

% For each kernel, compute prior
p_u = zeros(K, 1);

for k = 1:K
    k_win = max(1, floor(k-B/2)):min(K, floor(k+B/2));
    p_u(k) = sum((sqrt((u_sol(k) - u_sol(k_win)).^2 + 1e-10) ).^p);
end

% Scale according to lambda parameter
p_u = - p_u / (p * lambda^p);

end

