function f = eval_posterior( ...
                img_param, met_param, rf_data, alg, u_sol)
%EVAL_POSTERIOR
% Evaluates posterior probaility for a given data, prior and likelihood
% function. First, the total sample depth is separated in kernels of a
% given width and hop.
% size(rf_data) = K, N, M

% Evaluate likelihood function
like_val = eval_likelihood(img_param, rf_data, alg, u_sol);

% Evaluate prior
prior_val = eval_prior(img_param, met_param, u_sol);

% Accumulate log-probability
f = sum(like_val + prior_val);

end

