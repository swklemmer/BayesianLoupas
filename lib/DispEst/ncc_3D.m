function u = ncc_3D(est_p, RF_kern)
%NCC_3D Calculate initial displacement estimation for a collection of
% 3D kernels using the normalized cross-correlation. At every time step,
% the reference frame is the first frame in the ensemble. Output is
% expressed in samples.

% Retrieve signal spatial dimensions
len_z = size(RF_kern, 4);
len_x = size(RF_kern, 5);

% Preallocate displacement data
u = zeros(size(RF_kern, [1, 2, 3]));

% Generate 2D Tukey window
if est_p.crs_win
    tuck_win = tukeywin(len_z, .1) * tukeywin(len_x, .1)';
else
    tuck_win = 1;
end

for z = 1:size(RF_kern, 1)
    for x = 1:size(RF_kern, 2)
        for t = 1:size(RF_kern, 3)

            % Reference time-point
            R = 1 + est_p.crs_mref * (t - 1);

            % Get pre-deformation magnitude
            pre_kern = reshape(RF_kern(z, x, R, :, :, 1), [len_z, len_x]);

            % Get post-deformation magnitude
            post_kern = reshape(RF_kern(z, x, t, :, :, 2), [len_z, len_x]);
            
            % Calculate correlation (size = 2 w_z - 1 , 2 w_x - 1)
            win_corr = normxcorr2(tuck_win .* pre_kern, ...
                                  tuck_win .* post_kern);

            % Find maximum with coarse precission
            [~, max_corr] = max(win_corr, [], 'all');
            [max_z, ~] = ind2sub(size(win_corr), max_corr);

            % Calculate displacement
            u(z, x, t) = max_z - len_z;
        end
    end
end

% Optionally apply median filter
if est_p.crs_med
    u = reshape(medfilt2(reshape(real(u), size(u, 1), ...
        prod(size(u, 2:3))), [5 5]), size(u)); % [samples]
end

% If reference is a moving frame, accumulate over time
if est_p.crs_mref
    u = cumsum(u, 3); % [samples]
end

end
