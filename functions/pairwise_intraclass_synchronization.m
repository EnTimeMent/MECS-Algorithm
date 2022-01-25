function [Q, CM, C12, C21] = pairwise_intraclass_synchronization(time_series_1, time_series_2, class, tau, method, coincidence_function)

% Initializing outputs
CM = zeros(1, 1);
C12 = 0;
C21 = 0;
Q = 0;

% Checking input validity
if size(time_series_1, 1) > 1
    fprintf("Error: time series 1 should be a row vector\n");
    return
end
if size(time_series_2, 1) > 1
    fprintf("Error: time series 2 should be a row vector\n");
    return
end
if tau <= 0
    fprintf("Error: tau should be given as an integer number of samples\n");
    return
end

% Enumerating and counting the events of the given class in time series 1 and 2
event_indices_1 = find(time_series_1 == class);
event_indices_2 = find(time_series_2 == class);
m1 = size(event_indices_1, 2);
m2 = size(event_indices_2, 2);
if m1 == 0
    fprintf("Warning: time series 1 does not contain events of class %d\n", class);
    return
end
if m2 == 0
    fprintf("Warning: time series 2 does not contain events of class %d\n", class);
    return
end

% Computing coincidence matrix
CM = zeros(m1, m2);
for h = 1 : 1 : m1
    for k = 1 : 1 : m2
        x = event_indices_1(h);
        y = event_indices_2(k);
        switch coincidence_function
            case "constant"
                CM(h, k) = temporal_coincidence_constant(x, y, tau);
            case "linear"
                CM(h, k) = temporal_coincidence_linear(x, y, tau);
            otherwise
                CM(h, k) = temporal_coincidence_linear(x, y, tau);
        end
        
    end
end

% Computing temporal coincidence between time-series 1 and time-series 2
if method == "max"
    coincidence_scores_12 = max(CM, [], 2);
else
    coincidence_scores_12 = sum(CM, 2) ./ sum(CM ~= 0, 2);
end
coincidence_scores_12(isnan(coincidence_scores_12)) = 0;
C12 = sum(coincidence_scores_12) / m1;

% Computing temporal coincidence between time-series 2 and time-series 1
if method == "max"
    coincidence_scores_21 = max(CM, [], 1);
else
    coincidence_scores_21 = sum(CM, 1) ./ sum(CM ~= 0, 1);
end
coincidence_scores_21(isnan(coincidence_scores_21)) = 0;
C21 = sum(coincidence_scores_21) / m2;

% Computing intraclass pairwise synchronization
Q = (C12 + C21) / 2;

end