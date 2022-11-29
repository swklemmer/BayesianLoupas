function [Receive, TGC] = BMS_Receive(P, Trans)
% Specify structure arrays for reception.

% Specify Receive structure arrays
Receive = repmat(struct( ...
        'Apod', ones(1, Trans.numelements), ... % apodization
        'startDepth', P.startDepth, ...         % start depth [wvls]
        'endDepth', P.maxAcqLength,...          % end depth [wvls]
        'TGC', 1, ...                           % TCG structure
        'bufnum', 1, ...                        % RcvBuffer struct
        'framenum', 1, ...                      % RcvBuffer frame
        'acqNum', 0, ...                        % RcvBuffer adq
        'sampleMode', 'NS200BW', ...            % sample rate
        'mode', 0, ...                          % replace (0) or sum (1)
        'callMediaFunc', 1), ...                % enable scatterer movement
        1, P.bmode_adq);

% Assign a different adquisition to each structure
for i = 1:P.bmode_adq
    Receive(i).acqNum = i;
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [450, 590, 700, 800, 875, 925, 950, 975]; % Gain (0-1023)
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

end
