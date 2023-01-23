
% Iterate through each adquisition
adq_samples = size(RcvData{1}, 1)/P.bmode_adq - 128;

amplitude_data = zeros(adq_samples, 128);

for i = 1:P.bmode_adq
    RF_data = RcvData{1}(1 + (i - 1) * adq_samples : i * adq_samples, :);
    
    % Iterate over each column
    for j = 1:128
        break
        % Demodulate
        amplitude_data(:, j) = abs(hilbert(RF_data(:, j)));
    end

    amplitude_data = RF_data;

    imagesc(amplitude_data)
    pause()
end

close