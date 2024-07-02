function E = Eg(N, z, S0, S0z, S0zz, q0z, c, B, mu)
    % Initialize matrix
    E = zeros(2*N+1, 2*N+1) + 1i*zeros(2*N+1, 2*N+1);
    
    % Constants
    alpha = 1 ./ (1 + S0z.^2);
    f = S0z .* (q0z - c) .* alpha;
    gamma = c - q0z + f .* S0z;

    for mm = 1:(2*N+1)
        for jj = 1:(2*N+1)
            % Working with both code and matrix indices
            m = mm - 1;
            j = jj - 1;

            % Define integrand specific to matrix
            term = (gamma .* f + 3 .* S0zz .* S0z .* (alpha.^(5/2)) ...
                - ((alpha.^(3/2)) ./ S0) .* S0z - (alpha.^(3/2)) .* 1i .* (mu + m)) ...
                .* (1i .* (mu + m)) - (alpha.^(1/2)) ./ S0.^2 + B ./ S0.^3;
            
            trap = trapz(z, term .* exp(-1i * (m-j) * z));
            
            % Populate matrix entries
            E(jj, mm) = 1/(2*pi) * trap(1);
        end
    end
end

