function [p_xu, elapsed_t] = likelihood( ...
           img_param, met_param, alg, rf_lines, varargin)
%LIKELIHOOD
% Compute likelihood function using different methods:
% - IPS: Loupas likelihood with Interpolated Power Spectrum
% - ACK: Loupas likelihood with AutoCorrelation Kernel
% - NCC: Normalized Cross Correlation likelihood with spline interpolation

% Retrieve imaging parameters
f_c = img_param.f_c;
t_s = img_param.t_s;

% Retrieve method parameters
u_dim = met_param.u_dim;
ack_a = met_param.ack_a;
ncc_a = met_param.ncc_a;
qck_a = met_param.qck_a;

% Select algorithm
if strcmp(alg, 'ack')
    [p_xu, elapsed_t] = likelihood_ACK(...
                f_c, t_s, rf_lines, u_dim, ack_a, varargin);

elseif strcmp(alg, 'ncc')
    [p_xu, elapsed_t] = likelihood_NCC(...
                f_c, t_s, rf_lines, u_dim, ncc_a, varargin);

elseif strcmp(alg, 'qck')
    [p_xu, elapsed_t] = likelihood_QCK(...
                f_c, t_s, rf_lines, u_dim, qck_a, varargin);

elseif strcmp(alg, 'sqck')
    [p_xu, elapsed_t] = likelihood_SQCK(...
                f_c, t_s, rf_lines, u_dim, qck_a, varargin);
end

end
