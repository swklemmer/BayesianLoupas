function [img_param, met_param] = update_parameter(exp_param, value)
%UPDATE_PARAMETER

% Retrieve parameters
img_param = evalin('base', 'img_param');
met_param = evalin('base', 'met_param');

% Update corresponding parameter
if strcmp(exp_param, 'z_max')
    img_param.z_max = value;
    img_param.M = ceil(img_param.z_max / (img_param.f_c * img_param.t_s));

elseif strcmp(exp_param, 'SNR')
    img_param.SNR = value;

elseif strcmp(exp_param, 'N')
    img_param.N = value;

elseif strcmp(exp_param(1:5), 'alpha')
    met_param.ack_a = value;
    met_param.ncc_a = value;
end

end
