function A = Ag(N, z, S0z, q0z, c)
    % Initialize matrix
    A = zeros(2*N+1, 2*N+1) + 1i*zeros(2*N+1, 2*N+1);
    
    % Constants
    alpha = 1 ./ (1 + S0z.^2);
    f = S0z .* (q0z - c) .* alpha;

    % Loop over index to populate matrix
    for mm = 1:(2*N+1)
        for jj = 1:(2*N+1)
            % Working with both code and matrix indices
            m = mm - 1;
            j = jj - 1;

            % Define integrand specific to matrix
            term = f;

            trap = trapz(z, term .* exp(-1i * j * z));
            
            % Populate matrix entries
            A(jj, mm) = 1/(2*pi) * trap(1);
        end
    end
end