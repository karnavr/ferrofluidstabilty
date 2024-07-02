function B = Bg(N, z)
    % Initialize matrix
    B = zeros(2*N+1, 2*N+1) + 1i*zeros(2*N+1, 2*N+1);
    
    % Loop over index to populate matrix
    for mm = 1:(2*N+1)
        for jj = 1:(2*N+1)
            % Working with both code and matrix indices
            m = mm - 1;
            j = jj - 1;

            % Define integrand specific to matrix
            term = -1;
            
            trap = trapz(z, term .* exp(-1i * (m-j) * z));

            % Populate matrix entries
            B(jj, mm) = 1/(2*pi) * trap(1);
        end
    end
end