function C = Cg(N, z, S0, b, mu)
    % Initialize matrix
    C = zeros(2*N+1, 2*N+1) + 1i*zeros(2*N+1, 2*N+1);
    
    for mm = 1:(2*N+1)
        for jj = 1:(2*N+1)
            % Working with both code and matrix indices
            m = mm - 1;
            j = jj - 1;

            k = j + mu;

            % Define betas
            beta0 = beta(0, k, b, S0);
            beta1 = beta(1, k, b, S0);

            % Define integrand specific to matrix
            term = - (j + mu) .* S0 .* beta0;
            term = term ./ max(term);

            trap = trapz(z, term .* exp(-1i * (m - j) * z));

            % Populate matrix entries
            C(jj, mm) = 1/(2*pi) * trap(1);
        end
    end
end
