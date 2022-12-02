function compare_estimates(est_x, est_z, u_0, u_hat)
%COMPARE_ESTIMATES Summary of this function goes here

fig = figure();
subplot(1, 2, 1)
img1 = imagesc(est_x, est_z, squeeze(u_0(: ,:, 1)), ...
    [min(u_0, [], 'all'), max(u_0, [], 'all')]);
ylabel('Axial distance [\lambda]')
title('Initial estimation')
colorbar;

subplot(1, 2, 2)
img2 = imagesc(est_x, est_z, squeeze(u_hat(: ,:, 1)), ...
    [min(u_hat, [], 'all'), max(u_hat, [], 'all')]);
ylabel('Axial distance [\lambda]')
title('Bayesian estimation')
colorbar;

% Display all frames within 3 seconds
t = 1;
n = 0;
while n < 6
    pause(3 / size(u_0, 3))
    set(img1, 'CData', squeeze(u_0(:, :, t)));
    set(img2, 'CData', squeeze(u_hat(:, :, t)));

    if t < size(u_0, 3); t = t + 1;
        else; t = 1; n = n + 1; end
end

close(fig)

end

