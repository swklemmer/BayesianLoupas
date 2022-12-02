function eval_metric = likelihood_quality(p_x_u, u_dim, u_true, epsilon)
%LIKELIHOOD_QUALITY
% Integrates over small interval around true displacement

% Create integral dimension
[~, int_start] = find(u_dim >= u_true - epsilon, 1, 'first');
[~, int_end] = find(u_dim >= u_true + epsilon, 1, 'first');

% Calculate eval metric
eval_metric = sum(p_x_u(int_start:int_end));

end
