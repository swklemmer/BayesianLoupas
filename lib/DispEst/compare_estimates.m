function compare_estimates(est_x, est_z, u_1, u_2, u_3)
%COMPARE_ESTIMATES

% Minimum and maximum values
min_1 = min(u_1, [], 'all');
max_1 = max(u_1, [], 'all');
min_2 = min(u_2, [], 'all');
max_2 = max(u_2, [], 'all');
min_3 = min(u_3, [], 'all');
max_3 = max(u_3, [], 'all');

min_a = min([min_1, min_2, min_3]);
max_a = max([max_1, max_2, max_3]);

% Plot figures

fig = figure();
subplot(3, 2, 1)
img1 = imagesc(est_x, est_z, squeeze(u_1(: ,:, 1)), [min_a, max_a]);
xlabel('Lateral distance [\lambda]')
title('DAS Beamforming')
colorbar;

subplot(3, 2, 2)
img2 = imagesc(est_x, est_z, squeeze(u_2(: ,:, 1)), [min_a, max_a]);
ylabel('Axial distance [\lambda]')
title('FDAS-MAS Beamforming')
colorbar;

subplot(3, 2, 3)
img3 = imagesc(est_x, est_z, squeeze(u_3(: ,:, 1)), [min_a, max_a]);
title('True displacement')
colorbar;

subplot(3, 2, 4)
img4 = imagesc(est_x, est_z, squeeze(u_3(: ,:, 1)), [min_a, max_a]);
title('True displacement')
colorbar;

subplot(3, 2, 5)
img5 = imagesc(est_x, est_z, zeros(size(u_1, 1, 2)), [min_a, max_a]);
title('DAS Error')
colorbar;

subplot(3, 2, 6)
img6 = imagesc(est_x, est_z, zeros(size(u_1, 1, 2)), [min_a, max_a]);
title('FDAS-MAS Error')
colorbar;

% Display all frames within 10 seconds
t = 1;
n = 0;
while n < 6
    pause(10 / size(u_1, 3))

    set(img1, 'CData', squeeze(u_1(:, :, t)));
    set(img2, 'CData', squeeze(u_2(:, :, t)));
    set(img3, 'CData', squeeze(u_3(:, :, t)));
    set(img4, 'CData', squeeze(u_3(:, :, t)));
    set(img5, 'CData', squeeze(abs(u_1(:, :, t) - u_3(:, :, t))));
    set(img6, 'CData', squeeze(abs(u_2(:, :, t) - u_3(:, :, t))));

    if t < size(u_1, 3); t = t + 1;
        else; t = 1; n = n + 1; end
end

close(fig)

end

