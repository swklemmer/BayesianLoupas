function BFData = beamform_iq(IData, QData, fr_dec)

% Retrieve process parameters
bmode_adq = size(IData, 4);
n_ang = 3;

% Pre-allocate movie data
n_fr = floor(bmode_adq / fr_dec) - n_ang + 1; % frame number
BFData = zeros([size(IData, [1, 2]), n_fr]);

% Obtain magnitude of each frame
s_das = sqrt(squeeze(IData(:, :, 1, 1:fr_dec:end).^2 + ...
                     QData(:, :, 1, 1:fr_dec:end).^2));

for i=1:n_fr

    % Perform beamforming
    s_mas = zeros(size(IData, [1 2]));
    
    for m = i:(i+n_ang-2)
        for n = (m+1):(i+n_ang-1)
            s_mas = s_mas + sign(s_das(:, :, m) .* s_das(:, :, n)) .* ...
                sqrt(abs(s_das(:, :, m) .* s_das(:, :, n)));
        end
    end

    % Save envelope to ImageBuffer
    BFData(:, :, i) = s_mas;
end

end
