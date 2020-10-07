function max_mag = find_largest_magnitude(vec)

[maxV, idx] = max(abs(vec(:)));
max_mag = maxV * sign(vec(idx));