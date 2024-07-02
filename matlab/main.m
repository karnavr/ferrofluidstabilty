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

lambda = solveGenEig(solution, Nmodes, Nmu);  %% MAIN COMPUTATION

elapsedTime = toc; % stop timer and calculate elapsed time
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

% save the results to a file (name includes the number of modes and mu values)
save(['results/lambda_' num2str(Nmodes) '_' num2str(Nmu) '.mat'], 'lambda');

%% Plotting (can provide lambda as an argument to the function, or load it from a file)
% lambda = load('results/lambda_12_500.mat');
stabilityPlots(lambda);
