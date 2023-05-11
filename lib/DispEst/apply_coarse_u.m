function [I_disp, Q_disp] = apply_coarse_u(I_kern, Q_kern, u_coarse)
%APPLY_COARSE_U Displace the RF data by a discrete ammount of samples equal
% to the estimated coarse displacement. This way, the fine displacement
% estimator is responsible uniquely for the sub-sample displacement.
% INPUT: I_kern (size: [# Kernels in z, x, t], [Kernel size in z, x, t])
%        freqUS  (Center frequency)
%        FR      (Sampling frequency)
% OUTPUT: demod_data (size: z_len * x_len * t_len)

% Retrieve signal spatial dimensions
len_z = size(I_kern, 4);
len_x = size(I_kern, 5);

% Generate 2D Tukey window
tuck_win = tukeywin(len_z, .1) * tukeywin(len_x, .1)';

% Pre-allocate existing data
I_disp = I_kern;
Q_disp = Q_kern;

% Calculate differential disp
u_diff = cat(3, u_coarse(:, :, 1), diff(u_coarse, 1, 3));

for z = 1:size(I_kern, 1)
    for x = 1:size(I_kern, 2)
        for t = 1:size(I_kern, 3)

            % Apply Hanning window to second kernel
            I_kern_win = tuck_win .* ...
                reshape(I_kern(z, x, t, :, :, 2), [len_z, len_x]);
            Q_kern_win = tuck_win .* ...
                reshape(Q_kern(z, x, t, :, :, 2), [len_z, len_x]);

            % Displace second kernel circularly by the coarse estimation
            I_disp(z, x, t, :, :, 2) = circshift(I_kern_win, ...
                [real(u_diff(z, x, t)), imag(u_diff(z, x, t))]);
            Q_disp(z, x, t, :, :, 2) = circshift(Q_kern_win, ...
                [real(u_diff(z, x, t)), imag(u_diff(z, x, t))]);
        end
    end
end
end

