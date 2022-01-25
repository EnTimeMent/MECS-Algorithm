function [Q, S] = intraclass_synchronization(time_series, class, tau, method, coincidence_function)

% Checking input validity
if size(time_series, 2) > 1
    fprintf("Error: time series should be provided as a cell array having just one column\n");
    return
end
if size(time_series, 1) < 2
    fprintf("Error: at least two time series should be provided\n");
    return
end
if tau <= 0
    fprintf("Error: tau should be given as an integer number of samples\n");
    return
end

% Preparing data structures
num_time_series = size(time_series, 1);
P = nchoosek(1 : 1 : num_time_series, 2);
S = ones(num_time_series, num_time_series);
synchronization_values = zeros(size(P, 1), 1);

% Computing synchronization matrix
for k = 1 : 1 : size(P, 1)
    S(P(k, 1), P(k, 2)) = pairwise_intraclass_synchronization(...
        time_series{P(k, 1), 1}, time_series{P(k, 2), 1}, class, tau, method, coincidence_function);
    S(P(k, 2), P(k, 1)) = S(P(k, 1), P(k, 2));
    synchronization_values(k, 1) = S(P(k, 1), P(k, 2));
end

% Computing intraclass synchronization
Q = mean(synchronization_values);

end