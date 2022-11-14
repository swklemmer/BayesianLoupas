function p_xu = eval_likelihood(...
                        img_param, rf_data, alg, u_sol)
%EVAL_LIKELIHOOD
% EvAluate likelihood function at a given candidate solution u_sol using
% different methods:
% - ACK: Loupas likelihood with AutoCorrelation Kernel
% - NCC: Normalized Cross Correlation likelihood with spline interpolation
% - SSD: Sum of Squared Differences

% Select algorithm
if strcmp(alg, 'ack')
    p_xu = eval_likelihood_ACK(img_param, rf_data, u_sol);

elseif strcmp(alg, 'ncc')
    p_xu = eval_likelihood_NCC(img_param, rf_data, u_sol);

elseif strcmp(alg, 'ssd')
    p_xu = eval_likelihood_SSD(img_param, rf_data, u_sol);

end

end
