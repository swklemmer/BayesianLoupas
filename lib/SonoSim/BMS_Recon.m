function [Recon, ReconInfo] = BMS_Recon(P)
% Specify structure arrays for reconstruction.

% #1: IQ reconstruction into Interbuffer for displacement processing 
Recon(1) = struct( ...
    'senscutoff', 0.125, ...    % element angular sensitivity
    'pdatanum', 1, ...          % PixelData structure
    'rcvBufFrame', 1, ...       % RcvBuffer source frame
    'IntBufDest', [1, 1], ...   % Interbuffer dest. [buf, frame]
    'ImgBufDest', [0, 0], ...   % ImageBuffer dest. [buf, frame]
    'RINums', 1:P.bmode_adq);   % ReconInfo structure

    % Define ReconInfo structures.
    ReconInfo(1:P.bmode_adq) = repmat(struct( ...
        'mode', 'replaceIQ', ...        % data format
        'txnum', 0, ...                 % TX structure
        'rcvnum', 0, ...                % Receive structure
        'pagenum', 0, ...               % Destination page
        'regionnum', 1), ...            % PData region
        1, P.bmode_adq);

    % Assign a ReconInfo structure to each adquisition
    for i = 1:P.bmode_adq
        ReconInfo(i).txnum = 1 + mod(i, P.n_ang);
        ReconInfo(i).rcvnum = i;
        ReconInfo(i).pagenum = i;
    end

end
