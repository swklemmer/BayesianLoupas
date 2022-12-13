function f = eval_posterior_2D(met_param, rf_data, alg, u_sol)
%EVAL_POSTERIOR_3D
% Evaluates posterior probability for a given data, prior and likelihood
% function.
% The input data is passed as a 5D matrix:
% size(rf_data) = [# kernels in z, x; # samples in z, x & t]

% Evaluate likelihood function
like_val = eval_likelihood_2D(rf_data, alg, u_sol);

% Evaluate prior
prior_val = eval_prior_2D(met_param, u_sol);

% Accumulate log-probability
f = sum(like_val + prior_val, 'all');

end
