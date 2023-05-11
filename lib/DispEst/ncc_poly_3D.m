function u = ncc_poly_3D(est_p, RF_kern)
%NCC_3D Calculate initial displacement estimation for a collection of
% 3D kernels using the normalized cross-correlation with polynomial
% interpolation of sub-sample displacements.

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

% Create coarse and fine grids
[coarse_x, coarse_z] = meshgrid(-2:2, -2:2);
fine_dim = (-2:est_p.crs_sub:2);
[fine_x, fine_z] = meshgrid(fine_dim, fine_dim);

for z = 1:size(RF_kern, 1)
    for x = 1:size(RF_kern, 2)
        for t = 1:size(RF_kern, 3)

            % Reference time-point
            R = 1 + est_p.crs_mref * (t - 1);

            % Get pre-deformation kernel
            pre_kern = reshape(RF_kern(z, x, R, :, :, 1), [len_z, len_x]);

            % Get post-deformation kernel
            post_kern = reshape(RF_kern(z, x, t, :, :, 2), [len_z, len_x]);
            
            % Calculate correlation (size = 2 w_z - 1 , 2 w_x - 1)
            win_corr = normxcorr2(tuck_win .* pre_kern, ...
                                  tuck_win .* post_kern);

            % Find maximum with coarse precission
            [~, max_corr] = max(win_corr, [], 'all');
            [max_z, max_x] = ind2sub(size(win_corr), max_corr);

            % Establish coarse displacement limit (3, 2W - 3)
            max_z = min(max(max_z, 3), 2 * len_z - 3);
            max_x = min(max(max_x, 3), 2 * len_x - 3);

            % If kernel posseses less than 3 parrallel lines, use 1D fit
            if len_x < 3 
                % Fit correlation to polynomial around coarse maximum                
                xcorr_p = polyfit(-2:2, win_corr((-2:2) + max_z, 1), 4);

                % Evaluate polynomial with fine precission
                fine_xcorr = polyval(xcorr_p, fine_dim);

            else
                % Fit correlation to polynomial around coarse maximum
                xcorr_p = polyFit2D(...
                    win_corr((-2:2) + max_z, (-2:2) + max_x), ...
                    coarse_z, coarse_x, 4, 4);
    
                % Evaluate polynomial with fine precission
                fine_xcorr = polyVal2D(xcorr_p, fine_z, fine_x, 4, 4);
            end

            % Find maximum in fine grid
            [~, max_poly] = max(fine_xcorr, [], 'all');
            [~, max_dz] = ind2sub(size(fine_xcorr), max_poly);
  
            % Calculate displacement
            u(z, x, t) = max_z - len_z + fine_dim(max_dz);
        end
    end
end

end
