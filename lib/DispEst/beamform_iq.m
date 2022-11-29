function BFData = beamform_iq(IData, QData, fr_dec)

% Pre-allocate movie data
BFData = zeros(size(IData, [1, 2, 4]) ./ [1, 1, fr_dec] - [0, 0, -2]);
BFData = evalin('base', 'BFData;');

% Retrieve process parameters
bmode_adq = evalin('base', 'P.bmode_adq');
n_ang = evalin('base', 'P.n_ang');

% Obtain magnitude of each frame
s_das = sqrt(squeeze(IData(:, :, 1, :).^2 + QData(:, :, 1, :).^2));

for i=1:(bmode_adq-n_ang+1)

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

% Save data to workspace
assignin('base', 'BFData', BFData);

% Save data to workspace
assignin('base', 'IData', IData);
assignin('base', 'QData', QData);

close all

end
