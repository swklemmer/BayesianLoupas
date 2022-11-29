function Parameters = BMS_Parameters(P, Trans)
% Define system parameters.

Parameters = struct(...
    'numTransmit',      128, ...    % number of transmit channels.
    'numRcvChannels',   128, ...    % number of receive channels.
    'speedOfSound',     1540, ...   % speed of sound in m/sec
    'verbose',          0, ...
    'initializeOnly',   0, ...
    'simulateMode', P.simulate, ... % 0 = OFF, 1 = ON
    'waitForProcessing', 1, ...     % enable sync at end of Transfer
    'lambda', 0);                   % wavelength (to be calculated)

% Compute wavelength
Parameters.lambda = Parameters.speedOfSound / (Trans.frequency * 1e6);

end
