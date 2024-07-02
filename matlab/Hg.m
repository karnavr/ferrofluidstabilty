function H = Hg(N, z, S0, b, mu)
    % Initialize matrix
    H = zeros(2*N+1, 2*N+1) + 1i*zeros(2*N+1, 2*N+1);
    
    for mm = 1:(2*N+1)
        % Converting from code to matrix index
        m = mm;

        for jj = 1:(2*N+1)
            j = jj;

            k = mu + j;

            % Define beta
            beta0 = beta(0, k, b, S0);
            beta1 = beta(1, k, b, S0);

            % Define integrand specific to matrix
            term = S0 .* k .* (mu + m) .* beta1;

            trap = round(trapz(z, term .* exp(-1i * (m - j) * z)));
            
            % Populate matrix entries
            H(jj, mm) = 1/(2*pi) * trap(1);
        end
    end
end