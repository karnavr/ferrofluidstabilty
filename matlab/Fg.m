function F = Fg(N, z, S0, S0z, q0z, c, mu)
    % Initialize matrix
    F = zeros(2*N+1, 2*N+1) + 1i*zeros(2*N+1, 2*N+1);
    
    % Constants
    alpha = 1 ./ (1 + S0z.^2);
    f = S0z .* (q0z - c) .* alpha;
    gamma = c - q0z + f .* S0z;

    % Loop over index to populate matrix
    for mm = 1:(2*N+1)
        for jj = 1:(2*N+1)
            % Working with both code and matrix indices
            m = mm - 1;
            j = jj - 1;

            % Define integrand specific to matrix
            term = - gamma .* 1i .* (mu + m);

            trap = trapz(z, term .* exp(-1i * (m-j) * z));
            
            % Populate matrix entries
            F(jj, mm) = 1/(2*pi) * trap(1);
        end
    end
end


