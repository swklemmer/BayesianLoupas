function plot_mco(I_kern, Q_kern, I_moco, Q_moco, u_corr, z, x, t)
%PLOT_MCO Plot the result of motion correction at the kernel in position
% (z, x, t).

figure(9)
plot(squeeze(abs(I_kern(z, x, t, :, 1, 1) + ...
            1i * Q_kern(z, x, t, :, 1, 1))))
hold on
plot(squeeze(abs(I_kern(z, x, t, :, 1, 2) + ...
            1i * Q_kern(z, x, t, :, 1, 2))))
plot(squeeze(abs(I_moco(z, x, t, :, 1, 2) + ...
            1i * Q_moco(z, x, t, :, 1, 2))), '--')
hold off
grid on
legend({'Reference', 'Original frame', 'Corrected frame'})
title(sprintf("Motion correction by %d samples", ...
    u_corr(z, x, t+1) - u_corr(z, x, t)))

end

