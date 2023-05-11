function max_fc = pre_comp_SQCK(img_p, u_0, RF_kern)
%PRE_COMP_SQCK Pre-computes the maximum fast frequency column for a
%collection of kernels

% Find peak along frequency path
search_f = img_p.f_c * img_p.t_s + (-.1:1/(size(RF_kern, 4) * 10):.1);
rf_path = zeros(length(search_f), 1);

% Create DTFT dimentions
tau_m = 0:(size(RF_kern, 4) - 1);
tau_n = 0:(size(RF_kern, 6) - 1);

% Find maximum for each kernel
max_fc = zeros(size(RF_kern, [1 2]));

for z = 1:size(RF_kern, 1)
    for x = 1:size(RF_kern, 2)
            
        % Retrieve a single kernel
        rf_lines = squeeze(sum(RF_kern(z, x, 1, :, :, :), 5));

        for f = 1:length(search_f)
            % Calculate exponential mask at frequency path
            phi = exp(-2j * pi * (tau_m' * search_f(f) ...
                                + tau_n * search_f(f) * -u_0(z, x)));
        
            % Evaluate abs(DTFT) at desired point
            rf_path(f) = abs(rf_lines(:)' * phi(:));
        end

        % Find maximum index
        [~, max_index] = max(rf_path);
        max_fc(z, x) = search_f(max_index);
    end
end

end
