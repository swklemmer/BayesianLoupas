function compare_frame(est_x, est_z, u_0, u_hat, u_frame)
%COMPARE_ESTIMATES

z_lim = [min([min(u_0, [], 'all'), ...
    min(u_hat, [], 'all'), ...
    min(u_frame, [], 'all')]), ...
    max([max(u_0, [], 'all'), ...
    max(u_hat, [], 'all'), ...
    max(u_frame, [], 'all')])];

fig = figure();
fig.Position = [0, 400, 1100, 300];
subplot(1, 3, 1)
img1 = imagesc(est_x, est_z, u_0, z_lim);
ylabel('Axial distance [\lambda]')
title('Initial estimation')
colorbar;

subplot(1, 3, 2)
img2 = imagesc(est_x, est_z, u_hat, z_lim);
xlabel('Lateral distance [\lambda]')
title('Bayesian estimation')
colorbar;

subplot(1, 3, 3)
img3 = imagesc(est_x, est_z, u_frame, z_lim);
title('True displacement')
colorbar;

end

