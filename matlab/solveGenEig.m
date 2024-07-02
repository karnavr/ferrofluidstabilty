function lambda = solveGenEig(solution, Nmodes, Nmu, largest)

    % Default value for 'largest' parameter (because MATLAB doesn't support keyword arguments ugh)
    if nargin < 4
        largest = false;
    end

    % Constants
    Bond = 1.5;
    b = 0.1;
    epsilon = 1 - Bond/2;
    N = Nmodes;

    % Unpack solution
    coeffs = solution(2:end);
    c = solution(1);

    % Create domain and convert to real space
    z = linspace(-pi, pi, 100);
    [S0, S0z, S0zz] = fourierSeries(coeffs, z, pi);

    % Commonly used constants
    S0sq = 1 + S0z.^2;
    kappa = - (S0zz ./ (S0sq.^(3/2))) + (1 ./ (S0 .* S0sq.^(1/2)));
    q0z = c + (1 ./ S0) .* sqrt(S0sq .* ((c^2 + 2 * epsilon - 2 * kappa) .* (S0.^2) + Bond));

    % Set up stability stuff
    mu = linspace(0.001, 1.0, Nmu);

    % Create matrices that stay constant (don't depend on mu)
    A = Ag(N, z, S0z, q0z, c);
    B = Bg(N, z);
    D = zeros(2*N+1, 2*N+1);

    if largest
        lambda = zeros(length(mu), 1);

        for i = 1:Nmu
            % Create matrices
            E = Eg(N, z, S0, S0z, S0zz, q0z, c, Bond, mu(i));
            F = Fg(N, z, S0, S0z, q0z, c, mu(i));
            C = Cg(N, z, S0, b, mu(i));
            G = Gg(N, z, S0, S0z, q0z, b, c, mu(i));
            H = Hg(N, z, S0, b, mu(i));

            lhs = [A B; C D];
            rhs = [E F; G H];

            % Solve problem
            [sol, ~] = eigs(rhs, lhs, 1, 'largestreal');

            lambda(i) = sol(1);

            % Print progress
            if mod(i, 100) == 0
                fprintf('%d out of %d solved\n', i, Nmu);
            end
        end
    else
        lambda = zeros(4*N+2, Nmu);

        parfor i = 1:Nmu
            % Create matrices
            E = Eg(N, z, S0, S0z, S0zz, q0z, c, Bond, mu(i));
            F = Fg(N, z, S0, S0z, q0z, c, mu(i));
            C = Cg(N, z, S0, b, mu(i));
            G = Gg(N, z, S0, S0z, q0z, b, c, mu(i));
            H = Hg(N, z, S0, b, mu(i));

            lhs = [A B; C D];
            rhs = [E F; G H];

            % Solve problem
%             sol = eig(rhs, lhs);
            sol = eig(rhs, lhs, "qz"); % 

            % Save solution
            lambda(:, i) = sol;

        end
    end

    % Save solution to csv file
    % csvwrite(sprintf('stabilitySolutions/%d.%d.stabSol.csv', Nmu, Nmodes), lambda);

    return
end