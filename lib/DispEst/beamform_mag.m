function MData_bf = beamform_mag(img_param, IData, QData)

% Retrieve process parameters
bmode_adq = size(IData, 4);
n_ang = 3;

% Pre-allocate movie data
n_fr = floor(bmode_adq / img_param.fr_d) - n_ang + 1; % frame number
MData_bf = zeros([size(IData, [1, 2]), n_fr]);

% Decimate data in time
I_das = squeeze(IData(:, :, 1, 1:img_param.fr_d:end));
Q_das = squeeze(QData(:, :, 1, 1:img_param.fr_d:end));

% Add gaussian noise
I_n = db2mag(-img_param.snr) * rms(I_das, 'all');
I_das = I_das + I_n * randn(size(I_das));
Q_n = db2mag(-img_param.snr) * rms(Q_das, 'all');
Q_das = Q_das + Q_n * randn(size(Q_das));

% Take magnitude
s_das = sqrt(I_das.^2 + Q_das.^2);

% Perform beamforming
for i=1:n_fr
    s_mas = zeros(size(IData, [1 2]));
    
    for m = i:(i+n_ang-2)
        for n = (m+1):(i+n_ang-1)
            s_mas = s_mas + sign(s_das(:, :, m) .* s_das(:, :, n)) .* ...
                sqrt(abs(s_das(:, :, m) .* s_das(:, :, n)));
        end
    end

    % Save envelope to ImageBuffer
    MData_bf(:, :, i) = s_mas;
end

end
