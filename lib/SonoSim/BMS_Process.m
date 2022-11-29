function Process = BMS_Process()
% Specify Process structure arrays.

% # 1 : Non-Linear Beamforming process
Process(1) = struct( ...
    'classname', 'External', ...
    'method', 'BMS_NLbeamf', ...
    'Parameters', 0);

Process(1).Parameters = { ...
        'srcbuffer', 'inter', ...   % source buffer type
        'srcbufnum', 1, ...         % source buffer number
        'srcframenum', 1, ...       % source frame (0 = all)
        'srcsectionnum', 0, ...     % source section (0 = all)
        'srcpagenum', 0, ...        % source page (0 = all)
        'dstbuffer', 'none', ...    % destination buffer type
        'dstbufnum', 0, ...         % destination buffer number
        'dstframenum', 0, ...       % destination frame
        'dstsectionnum', 0, ...     % destination section
        'dstpagenum', 0};           % destination page

end

% # 1: Display last B-mode
% Process(1) = struct( ...
%     'classname', 'Image', ...
%     'method', 'imageDisplay', ...
%     'Parameters', 0);
% 
% Process(1).Parameters = { ...
%         'imgbufnum', 1, ...            % source ImageBuffer
%         'framenum', 1, ...             % source frame
%         'pdatanum', 1, ...             % PixelData structure
%         'pgain', 1, ...                % digital gain
%         'reject', 0, ...               % reject level (0-100)
%         'persistMethod', 'simple', ... % persist mode
%         'persistLevel', 0, ...         % persist level (0-100)
%         'interpMethod', '4pt', ...     % interp. method
%         'grainRemoval', 'none', ...    % spatial filter
%         'processMethod', 'none', ...   % speckle reduction
%         'averageMethod', 'none', ...   % temporal filter
%         'compressMethod', 'log', ...   % compression mode 
%         'compressFactor', 15, ...      % comp. level (0-100)
%         'mappingMethod', 'full', ...   % gray scale range
%         'display', 1,...               % enable display
%         'displayWindow', 1};           % DispWindow structure