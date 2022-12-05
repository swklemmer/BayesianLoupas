function compare_frame(est_x, est_z, u_0, u_hat, u_frame)
%COMPARE_ESTIMATES

fig = figure();
subplot(1, 3, 1)
img1 = imagesc(est_x, est_z, u_0, ...
    [min(u_0, [], 'all'), max(u_0, [], 'all')]);
ylabel('Axial distance [\lambda]')
title('Initial estimation')
colorbar;

subplot(1, 3, 2)
img2 = imagesc(est_x, est_z, u_hat, ...
    [min(u_hat, [], 'all'), max(u_hat, [], 'all')]);
ylabel('Axial distance [\lambda]')
title('Bayesian estimation')
colorbar;

subplot(1, 3, 3)
img3 = imagesc(est_x, est_z, u_frame, ...
    [min(u_frame, [], 'all'), max(u_frame, [], 'all')]);
ylabel('Axial distance [\lambda]')
title('True displacement')
colorbar;

end

