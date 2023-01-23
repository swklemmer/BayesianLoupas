met_param = struct(...
    'u_dim', -0.5:1e-3:0.5, ... % [lmbds]
    'ack_a', 5e-3, ...           % [distribution width]
    'ncc_a', 2e-2);              % [distribution width]

likelihood_pdf()
figure(6)
hold on
met_param.ack_a = 1e-2;
likelihood_pdf
met_param.ack_a = 5e-2;
likelihood_pdf
met_param.ack_a = 1e-1;
likelihood_pdf
met_param.ack_a = 5e-1;
likelihood_pdf
hold off

figure(7)
close

figure(8)
close

fig = figure(6);
fig.Position = [1000, 200, 300, 200];
ylim([0, 0.07])
%xlim([-0.43 0.03])
xlim(u_true + [-0.05 0.05])
legend({'\beta = 5e-3', '', '\beta = 1e-2', '', '\beta = 5e-2', '', ...
    '\beta = 1e-1', '', '\beta = 5e-1', 'True'}, ...
    'FontSize', 10, 'Location', 'northwest')
title('ACK likelihood PDF', 'Fontsize', 12)
xlabel('Displacement [\lambda]', 'Fontsize', 12)
ylabel('Probability')
title('ACK likelihood PDF', 'Fontsize', 14)