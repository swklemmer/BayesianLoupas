function e_met = error_metrics(u_hat, u_tru, mean_u)
%ERROR_METRICS Calculates mean bias, variance and root square error of
%displacement estimation.

e_met = zeros(1, 3);

e_met(1) = mean(u_hat - u_tru, 'all') / mean_u;
e_met(2) = var(u_hat - u_tru, [], 'all') / mean_u;
e_met(3) = sqrt(e_met(2) + e_met(1)^2);

end
