function BMS_NLbeamf(IData, QData)
%BMS_NLBEAMF

% Declare variables common to all processes
persistent bmode_adq n_ang BFData

% Run only at initialization
if isempty(bmode_adq)

    % Retrieve process parameters
    bmode_adq = evalin('base', 'P.bmode_adq');
    n_ang = evalin('base', 'P.n_ang');

    % Pre-allocate movie data
    evalin('base', sprintf('BFData = zeros(%d, %d, %d);', ...
        size(IData, 1), size(IData, 2), bmode_adq-n_ang+1));
    BFData = evalin('base', 'BFData;');

    % Save estimation dimentions to workspace
    PData = evalin('base', 'PData');
    assignin('base', 'img_x', (0:PData.Size(2)-1) * PData.PDelta(1));
    assignin('base', 'img_z', (0:PData.Size(1)-1) * PData.PDelta(3));
end

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
