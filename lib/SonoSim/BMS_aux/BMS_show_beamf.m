function BMS_show_beamf(hObject, ~, varargin)
%BMS_SHOW_MOVIE
    persistent fig

    % Retrieve estimated displacement and dimentions
    BFData = evalin('base', 'BFData');
    img_x = evalin('base', 'img_x');
    img_z = evalin('base', 'img_z');

    fig = figure();
    img = imagesc(img_x, img_z, squeeze(BFData(: ,:, 1)), ...
        [0, 2e9]);
    colormap('jet')
    title(sprintf('Beamformed B-modes with C_t = %4.2f [m/s]', varargin{1}))
    xlabel('Lateral direction [\lambda]')
    ylabel('Axial depth [\lambda]')
    colorbar;

    % Display all frames within 3 seconds
    t = 1;
    while hObject.Value
        pause(3 / size(BFData, 3))

        set(img, 'CData', squeeze(BFData(:, :, t)));
        if t < size(BFData, 3); t = t + 1;
        else; t = 1; break; end
    end

    close(fig)
end
