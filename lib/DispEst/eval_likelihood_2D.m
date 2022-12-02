function p_xu = eval_likelihood_2D(rf_data, alg, u_sol)
%EVAL_LIKELIHOOD_3D
% Evaluate likelihood function at a given candidate solution u_sol using
% different methods:
% - ACK: Loupas likelihood with AutoCorrelation Kernel
% - NCC: Normalized Cross Correlation likelihood with spline interpolation
% - SSD: Sum of Squared Differences

% Select algorithm
if strcmp(alg, 'ack')
    p_xu = eval_likelihood_ACK_2D(rf_data, u_sol);

elseif strcmp(alg, 'ncc')
    p_xu = eval_likelihood_NCC_2D(rf_data, u_sol);

elseif strcmp(alg, 'ssd')
    p_xu = eval_likelihood_SSD_2D(rf_data, u_sol);

end

end
