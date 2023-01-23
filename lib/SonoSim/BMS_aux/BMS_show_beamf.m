function BMS_show_beamf(hObject, ~, varargin)
%BMS_SHOW_MOVIE
    persistent fig

    % Retrieve estimated displacement and dimentions
    BFData = evalin('base', 'BFData');
    img_x = evalin('base', 'img_x');
    img_z = evalin('base', 'img_z');

    fig = figure();
    fig.Position = [100, 100, 250, 250];
    img = imagesc(img_x, img_z, squeeze(BFData(: ,:, 1)), ...
        [min(BFData, [], 'all'), max(BFData, [], 'all')]);
    colormap('jet')
    if ~isempty(varargin)
        title(sprintf('Beamformed B-modes with C_t = %4.2f [m/s]', varargin{1}))
    else
        %title('Beamformed B-modes', 'Interpreter','latex')
    end
    xlabel('Lateral direction [$\lambda$]', 'Interpreter','latex')
    ylabel('Axial depth [$\lambda$]', 'Interpreter','latex')
    %colorbar;

    % Display all frames within 3 seconds
    t = 1;
    n = 0;
    while n < 1
        %pause(3 / size(BFData, 3))
    
        set(img, 'CData', squeeze(BFData(:, :, t)));
        if t < size(BFData, 3); t = t + 1;
        else; t = 1; n = n + 1; end

        pause()
    end

    close(fig)
end
