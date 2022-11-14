function [rf_pulse, t_cell] = create_pulse(img_param, varargin)
%CREATE_PULSE Returns time steps, complex IQ signal and -6dB resolution
%time of RF pulse.

% Retrieve imaging parameters
f_c = img_param.f_c;
bw = img_param.bw;
t_s = img_param.t_s;

% Create time steps
N_pulse = 2 * ceil(1e-6 / t_s / 2); % [smpls]
t_pulse = ((1:N_pulse) - N_pulse/2) * t_s; % [s]

% Create complex IQ signal
[i_pulse, q_pulse] = gauspuls(t_pulse, f_c, bw);
rf_pulse = complex(i_pulse, q_pulse);

% Calculate -6dB resolution cell
t_cell = gauspuls('cutoff', f_c, bw, -6, -6) * 2; % [s]

% Show pulse
if size(varargin) > 0
    fig = figure(1);
    fig.Position = [0, 600, 300, 200];
    plot(t_pulse * 1e6, abs(rf_pulse), '--');
    hold on
    plot(t_pulse * 1e6, real(rf_pulse), '--');
    grid on;
    xlabel('Time [us]')
    title('RF Pulse')
end

end
