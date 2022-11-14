addpath('lib/')

% Experiment info
alg_list = {'ack'};
exp_param = 'alpha';
param_list = logspace(-2.5, 0, 20);
N_list = 2:2:16;

% Imaging parameters
img_param = struct(...
    'f_c',      7e6, ...    % [Hz] 
    'bw',       0.5, ...    % [%]
    't_s',      0, ...      % [s]
    'z_max',    5, ...     % [lmbds]
    'M',        10, ...      % [number of samples]
    'N',        3, ...      % [ensemble length]
    'SNR',      10, ...     % [dB]
    'SNR_rho',  1);         % [peak to peak]

img_param.t_s = 1 / (4 * img_param.f_c);
img_param.M = ceil(img_param.z_max / (img_param.f_c * img_param.t_s));
img_param.SNR_rho = 10^(img_param.SNR/20);

% Method parameters
met_param = struct(...
    'u_dim', -0.5:0.001:0.5, ... % [lmbds]
    'ack_a', 0.01, ...           % [distribution width]
    'ncc_a', 1);              % [distribution width]

% Simulation parameters
rng(420)
n_exp = 2e1; % [number of experiments]
u_true = max(min(randn([1, n_exp]) / 3 , 0.5), -0.5); % [lmbds]
int_val = zeros(length(param_list), length(N_list), n_exp);

% Create gaussian RF Pulse 
[rf_pulse, t_cell] = create_pulse(img_param);

for e = 1:length(N_list)
    e
    % Assign value to corresponing parameter
    [img_param, met_param] = update_parameter('z_max', N_list(e));

    for p = 1:length(param_list)
        % Assign value to corresponing parameter
        [img_param, met_param] = update_parameter(exp_param, param_list(p));
    
        for alg = 1:length(alg_list)
    
            % Raise parameter change flag
            param_flag = 1;
    
            for n = 1:n_exp % [lmbds]
        
                % Create RF lines
                rf_lines = create_rf_line(...
                    img_param, rf_pulse, t_cell, u_true(n));
    
                % Calculate likelihood function using selected algorithm
                [p_xu, ~] = likelihood(...
                    img_param, met_param, alg_list{alg}, rf_lines);
    
                % Integrate probability over +- lambda/2
                int_val(p, e, n) = sum(p_xu) * (met_param.u_dim(2) - met_param.u_dim(1));
            end
        end
    end
end

%% Process results
int_mean = mean(int_val, 3);
%f_MN = sqrt(1 ./ abs(int_mean - 1)) / sqrt(met_param.ack_a);

% n2_mean = int_mean(: , 1);
% n3_mean = int_mean(: , 2);
% n4_mean = int_mean(: , 3);
% n5_mean = int_mean(: , 4);
% n6_mean = int_mean(: , 5);
% n7_mean = int_mean(: , 6);
% n8_mean = int_mean(: , 7);
% n9_mean = int_mean(: , 8);
% n10_mean = int_mean(: , 9);
% n11_mean = int_mean(: , 10);
% n12_mean = int_mean(: , 11);

int_std = std(int_val, [], 3);

%%
figure(1)
mesh(param_list, N_list, int_mean');
zlim([0, 2])
xlabel('\alpha')
ylabel('N')

h= gca;
set(h,'xscale','log')
%set(h,'zscale','log')
% 
% %%
% figure(2)
% mesh(param_list, N_list, int_std' ./ int_mean');
% zlim([-1, 1])
% xlabel('\alpha')
% ylabel('N')

%%
% figure(3)
% semilogx(param_list, int_mean)
% ylim([0, 2])
% legend(cellstr(num2str(N_list', 'Z_{max}=%d')))
%%
%save(sprintf('Z%d_N%d.mat', img_param.z_max, img_param.N), 'img_param', 'met_param', 'int_mean', 'int_std', 'n_exp')
