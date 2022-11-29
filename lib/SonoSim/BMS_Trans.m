function [P, Trans] = BMS_Trans(P)
% Define Trans structure array.

Trans = struct( ...
    'name',             'L11-5v', ...      % known transducer name
    'units',            'wavelengths', ... % distance units
    'maxHighVoltage',   90);               % pulser voltage limit

Trans = computeTrans(Trans);

% Calculate maximum adquisition length [wvls]
P.maxAcqLength = ... % longest distance between element and scatterer
    ceil(sqrt(P.endDepth^2 + ((Trans.numelements - 1) * Trans.spacing)^2));

end
