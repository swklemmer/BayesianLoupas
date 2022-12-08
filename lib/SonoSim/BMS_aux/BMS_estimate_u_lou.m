function BMS_estimate_u_lou(IData, QData)
%BMS_ESTIMATE_U_LOU: Loupas algorithm with median filtering.
% Input: IData, QData
% Size(IData) =
%   (Z dim, X dim, frame, page) = (PData.Size([1, 2]), 1, P.bmode_adq)
% Output: MovieData
% Size(MovieData) = (Z est, X est, P.bmode_adq-1) 
% 
% PData.Size([1,2]) = [140, 43] @ pdelta = 0.3

% Declare static variables
persistent bmode_adq MovieData ens_len med_sz w_x w_z hop_x hop_z N_x N_z ...
           hann_win_F hann_win_f cum_sum
    
% Read parameter change flag
param_flag = evalin('base', 'param_flag');

% Run only every time parameters change
if param_flag

    % Lower flag
    assignin('base', 'param_flag', 0);

    % Retrieve constant parameters
    bmode_adq = evalin('base', 'P.bmode_adq'); % Number of adquisitions

    % Get ROI dimensions
    x_max = evalin('base', 'PData.Size(2)'); % ROI size [wvls]
    z_max = evalin('base', 'PData.Size(1)'); % ROI size [wvls]
    dx = evalin('base', 'PData.PDelta(1)');  % z resolution [wvls]
    dz = evalin('base', 'PData.PDelta(3)');  % z resolution [wvls]

    % Get estimation Parameters (may change due to grid resolution)
    est_param = evalin('base', 'current_param.lou');
    ens_len = est_param.ens_len; % N = Ensemble length
    axi_len = est_param.axi_len; % Axial window length [wvls]
    axi_hop = est_param.axi_hop; % Axial window hop [wvls]
    lat_len = est_param.lat_len; % Lateral window length [wvls]
    lat_hop = est_param.lat_hop; % Lateral window hop [wvls]
    med_sz = est_param.med_sz;  % Median filter size [smpls]
    cum_sum = est_param.cum_sum; % Moving average length [smpls]
    
    % Define window parameters
    w_x = ceil(lat_len / dx);            % Lateral window length [smpls]
    w_z = ceil(axi_len / dz);            % M = Axial window length [smpls]
    hop_x = max(floor(lat_hop / dx), 1); % Lateral hop size
    hop_z = max(floor(axi_hop / dz), 1); % Axial hop size [smpls]
    N_x = floor((x_max - w_x) / hop_x);  % Number of lateral windows
    N_z = floor((z_max - w_z) / hop_z);  % Number of axial windows

    % Save parameters to text file
    param = {'Ensemble Length'; 'Axi. Win. Size'; 'Axi. Win. Hop'; ...
            'Lat. Win. Size'; 'Lat. Win. Hop'; 'Median Filter Size'};
    value = [ens_len; w_z * dz; hop_z * dz; w_x * dx; hop_x * dx; med_sz];
    units = {'smpls'; 'wvls'; 'wvls'; 'wvls'; 'wvls'; 'smpls'};
    param_table = table(param, value, units);
    evalin('base', 'param_table = table();');
    assignin('base', 'param_table', param_table);

    % Generate estimation dimentions [wvls]
    evalin('base', ...
        sprintf(...
        ['est_x = PData.Origin(1) + (%d/2 + (0:(%d-1)) * %d)', ...
        '* PData.PDelta(1);'], w_x, N_x, hop_x))

    evalin('base', ...
        sprintf(...
        ['est_z = PData.Origin(3) + (%d/2 + (0:(%d-1)) * %d)',...
        '* PData.PDelta(3);'], w_z, N_z, hop_z))

    % Create hann window
    %hann_win_F = repmat(hann(w_z), 1, w_x);
    %hann_win_f = repmat(hann(w_z-1), 1, w_x);
    hann_win_F = ones(w_z, w_x);
    hann_win_f = ones(w_z-1, w_x);

    % Preallocate displacement
    evalin('base', sprintf('MovieData = zeros(%d, %d, %d);', ...
                            N_z, N_x, bmode_adq-ens_len+1))
    MovieData = evalin('base', 'MovieData');
end

% Estimate velocity
for t = 1:(bmode_adq-ens_len+1)
    % Calculate ensemble windows
    win_t_F = (1:ens_len-1) + (t - 1);
    win_t_f = (1:ens_len) + (t - 1);

    for x = 1:N_x
        % Calculate lateral window
        win_x = (1:w_x) + (x - 1) * hop_x;

        for z = 1:N_z
            % Calculate axial windows
            win_z_F = (1:w_z) + (z - 1) * hop_z;
            win_z_f = (1:w_z-1) + (z - 1) * hop_z;

            % Calculate frecuency estimates
            F_0 = atan(...
                sum(hann_win_F .* ...
                (QData(win_z_F, win_x, 1, win_t_F) .* ...
                IData(win_z_F, win_x, 1, win_t_F + 1) - ...
                IData(win_z_F, win_x, 1, win_t_F) .* ...
                QData(win_z_F, win_x, 1, win_t_F + 1)), 'all') / ...
                sum(hann_win_F .* ...
                (IData(win_z_F, win_x, 1, win_t_F) .* ...
                IData(win_z_F, win_x, 1, win_t_F + 1) + ...
                QData(win_z_F, win_x, 1, win_t_F) .* ...
                QData(win_z_F, win_x, 1, win_t_F + 1)), 'all')) / (2*pi);
        
            f_dopp = atan(...
                sum(hann_win_f .* ...
                (QData(win_z_f, win_x, 1, win_t_f) .* ...
                IData(win_z_f + 1, win_x, 1, win_t_f) - ...
                IData(win_z_f, win_x, 1, win_t_f) .* ...
                QData(win_z_f + 1, win_x, 1, win_t_f)), 'all') / ...
                sum(hann_win_f .* ...
                (IData(win_z_f, win_x, 1, win_t_f) .* ...
                IData(win_z_f + 1, win_x, 1, win_t_f) + ...
                QData(win_z_f, win_x, 1, win_t_f) .* ...
                QData(win_z_f + 1, win_x, 1, win_t_f)), 'all')) / (2*pi);
        
            % Calculate differential displacement [wvls]
            MovieData(z, x, t) =  F_0 / (1 + f_dopp);
        end
    end

    % Apply median filter
    %MovieData(:, :, t) = medfilt2(MovieData(:, :, t), ...
    %                    [med_sz, med_sz], 'symmetric');
end

% Save accumulated displacement [wvls] to workspace 
assignin('base', 'MovieData', ...
    convn(MovieData, ones(1, 1, cum_sum)/cum_sum, 'same'));

end
