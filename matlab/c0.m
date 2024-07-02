function c0_val = c0(k, b, B)
    % Wave speed for small amplitude waves, depending on the wave-number k
    c0_val = sqrt((1 ./ k) .* ((-beta(1, k, b, 1) ./ beta(0, k, b, 1)) .* (k.^2 - 1 + B)));
end

