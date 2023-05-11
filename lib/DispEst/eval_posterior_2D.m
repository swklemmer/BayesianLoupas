function f = eval_posterior_2D(img_p, met_p, rf_data, u_sol)
%EVAL_POSTERIOR_3D
% Evaluates posterior probability for a given data, prior and likelihood
% function.
% The input data is passed as a 5D matrix:
% size(rf_data) = [# kernels in z, x; # samples in z, x & t]

% Evaluate likelihood function
like_val = eval_likelihood_2D(img_p, met_p.alg, rf_data, u_sol);

% Evaluate prior
prior_val = eval_prior_2D(met_p, u_sol);

% Accumulate log-probability
f = sum(like_val + prior_val, 'all');

% fprintf('Like = %.1f | Prior = %.1f | Prop. = %.1f\n', ...
%     sum(like_val, 'all'), sum(prior_val, 'all'), ...
%     100 * abs(sum(prior_val, 'all') / sum(like_val, 'all')))
end
