%% Solving the generalized eigenvalue problem

% import full solution branch 
solution_branch = load('32.1.0e-15.100.Broyden/solutions_32.1.0e-15.100.Broyden.dat');

% take a solution from somewhere in the middle
solution = solution_branch(24,:);

% set parameters
Nmodes = 24;       % number of perturbation modes
Nmu = 1000;        % number of mu values between 0 and 1

%% solve the problem (and time it)
tic; % start the timer

lambda = solveGenEig(solution, Nmodes, Nmu, true);  %% MAIN COMPUTATION

elapsedTime = toc; % stop timer and calculate elapsed time
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

% save the results to a file (name includes the number of modes and mu values)
save(['results/lambda_' num2str(Nmodes) '_' num2str(Nmu) '.mat'], 'lambda');

%% Plotting (can provide lambda as an argument to the function, or load it from a file)
% lambda = load('results/lambda_12_500.mat');
stabilityPlots(lambda);

%%

mu = 0.5;

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


A = Ag(N, z, S0z, q0z, c);
B = Bg(N, z);
D = zeros(2*N+1, 2*N+1);
E = Eg(N, z, S0, S0z, S0zz, q0z, c, Bond, 0.5);
F = Fg(N, z, S0, S0z, q0z, c, mu);
C = Cg(N, z, S0, b, mu);
G = Gg(N, z, S0, S0z, q0z, b, c, mu);
H = Hg(N, z, S0, b, mu);

% surf(real([A B; C D]))
% 
% figure
% surf(imag([A B; C D]))

% plot matrix with colorbar
figure
imagesc(real(C(1:end-5,:)))

%% 
test = beta(0, 50, 0.1, ones(1, 100));
plot(test, 'o')
axis tight

