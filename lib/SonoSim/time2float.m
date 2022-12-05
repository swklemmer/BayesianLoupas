function [sorted_times, sort_ind] = time2float(adq_times)
% Auxiliar Function for parsing adquisition times in h5 file

% Convert string to floats
float_times = zeros(1, length(adq_times));

for i = 1:length(adq_times)
    t_string = strrep(char(adq_times(i)), '_', '.');
    float_times(i) = str2double(t_string);
end

[sorted_times, sort_ind] = sort(float_times);
end

