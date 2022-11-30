function Data_bf = beamform_rf(img_param, TX, IData, RFData)



% Retrieve process parameters
bmode_adq = size(IData, 4);
smpls_fr = size(RFData, 1) / bmode_adq - 128; % samples per frame
n_ang = 3; % steering angle number
n_fr = floor(bmode_adq / img_param.fr_d) - n_ang + 1; % frame number

% Split Data into acquisitions
AData = zeros(smpls_fr, size(RFData, 2), bmode_adq);

for i = 1:bmode_adq
    % Delay every element according to transmit steering angle
    elem_t = round(TX(1 + mod(i, n_ang)).Delay * 4); % one wvl = 4 spls
    
    for j = 1:size(RFData, 2)
        AData(:, j, i) = ...
            RFData((1:smpls_fr) + (i - 1) * smpls_fr + elem_t(j), j);
    end
end

% Add gaussian noise
n_pow = db2mag(-img_param.snr) * rms(AData, 'all');
s_das = AData + n_pow * randn(size(AData));

% Design hanning windowed low-pass filter
taps_lp = fir1(501, 0.15, 'low');

% Perform beamforming
Data_bf = zeros([size(AData, [1, 2]), n_fr]);

for i=1:n_fr
    s_mas = zeros(size(AData, [1 2]));
    
    for m = i:(i+n_ang-2)
        for n = (m+1):(i+n_ang-1)
            s_mas = s_mas + sign(s_das(:, :, m) .* s_das(:, :, n)) .* ...
                sqrt(abs(s_das(:, :, m) .* s_das(:, :, n)));
        end
    end

    % Save filtered data to buffer
    Data_bf(:, :, i) = filtfilt(taps_lp, 1, s_mas);
end

end
