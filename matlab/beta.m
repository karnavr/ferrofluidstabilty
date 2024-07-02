function beta_val = beta(n, k, b, S0)
    beta1 = besseli(1, k*b) .* besselk(n, k.*S0);
    beta2 = (-1)^n .* besselk(1, k*b) .* besseli(n, k.*S0);

    beta_val = beta1 + beta2;
end