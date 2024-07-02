function G = Gg(N, z, S0, S0z, q0z, b, c, mu)
    % Initialize matrix
    G = zeros(2*N+1, 2*N+1) + 1i*zeros(2*N+1, 2*N+1);
    
    for mm = 1:(2*N+1)
        % Converting from code to matrix index
        m = mm;

        for jj = 1:(2*N+1)
            j = jj;

            k = mu + j;

            % Define betas
            beta0 = beta(0, k, b, S0);
            beta1 = beta(1, k, b, S0);

            % Define integrand specific to matrix
            term = S0 .* S0z .* c .* k.^2 .* beta1 ...
                 - S0z .* c .* k .* beta0 ...
                 + 1i .* S0 .* q0z .* k.^2 .* beta0 ...
                 - 1i .* S0 .* c .* k .* (mu + m) .* beta0;

            trap = trapz(z, term .* exp(-1i * (m - j) * z));
            
            % Populate matrix entries
            G(jj, mm) = 1/(2*pi) * trap(1);
        end
    end
end