function c = temporal_coincidence_constant(x, y, tau)
d = temporal_distance(x, y);
if d < tau
    c = 1;
else
    c = 0;
end