function c = temporal_coincidence_linear(x, y, tau)
d = temporal_distance(x, y);
if d < tau
    c = 1 - (d / tau);
else
    c = 0;
end