function UI = BMS_GUI(P, Parameters, DisplayWindow, Recon)
% Specify UI Control Elements

% - Sensitivity Cutoff
UI(1) = struct( ...
    'Control', 0, ...
    'Callback', {@SensCutoffCallback});

UI(1).Control = {'UserB7', ...
                'Style', 'VsSlider', ...
                'Label', 'Sens. Cutoff', ...
                'SliderMinMaxVal', [0, 1.0, Recon(1).senscutoff], ...
                'SliderStep', [0.025, 0.1], ...
                'ValueFormat','%1.3f'};

% - Range Change
MinMaxVal = [64, 300, P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(DisplayWindow(1),'AxesUnits')&&~isempty(DisplayWindow(1).AxesUnits)
    if strcmp(DisplayWindow(1).AxesUnits, 'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * Parameters.lambda * 1e3;
    end
end

UI(2) = struct( ...
    'Control', 0, ...
    'Callback', {@TDE_show_movie});

UI(2).Control = {'UserC1', ...
                 'Style', 'VsToggleButton', ...
                 'Label', 'Play SW loop'};

end

% **** Callback routines ****

% SensCutoffCallback - Sensitivity cutoff change
function SensCutoffCallback(hObject, ~)
    
    % Get slider value
    UIValue = get(hObject, 'Value');

    % Change sensitivity cutoff in Recon structure
    ReconL = evalin('base', 'Recon');
    for i = 1:size(ReconL, 2)
        ReconL(i).senscutoff = UIValue;
    end
    assignin('base', 'Recon', ReconL);

    % Change control Structure to force an update
    Control = evalin('base', 'Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon'};
    assignin('base', 'Control', Control);
    return
end
