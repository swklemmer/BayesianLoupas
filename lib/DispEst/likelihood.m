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
M = img_param.M;
N = img_param.N;

% Retrieve method parameters
u_dim = met_param.u_dim;
ack_a = met_param.ack_a;
ncc_a = met_param.ncc_a;

% Select algorithm
if strcmp(alg, 'ips')
    [p_xu, elapsed_t] = likelihood_IPS(...
                f_c, t_s, M, N, rf_lines, u_dim, varargin);

elseif strcmp(alg, 'ack')
    [p_xu, elapsed_t] = likelihood_ACK(...
                f_c, t_s, M, N, rf_lines, u_dim, ack_a, varargin);

elseif strcmp(alg, 'ncc')
    [p_xu, elapsed_t] = likelihood_NCC(...
                f_c, t_s, N, rf_lines, u_dim, ncc_a, varargin);

elseif strcmp(alg, 'ssd')
    [p_xu, elapsed_t] = likelihood_SSD(...
                f_c, t_s, M, N, rf_lines, u_dim, ncc_a, varargin);
end

end
