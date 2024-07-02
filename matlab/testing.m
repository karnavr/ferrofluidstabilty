%% Testing fourierSeries()

% Coefficients for cosine
coefficients = [0; 1]; 

% Define the domain
domain = linspace(-pi, pi, 100);
L = pi;

% Call the fourierSeries function
[S, Sz, Szz] = fourierSeries(coefficients, domain, L);

% Plot
figure;
subplot(3, 1, 1);
plot(domain, S, 'b', 'LineWidth', 2);
title('Fourier Series: S(x)');
xlabel('x');
ylabel('S(x)');
grid on;

subplot(3, 1, 2);
plot(domain, Sz, 'r', 'LineWidth', 2);
title('First Derivative: S''(x)');
xlabel('x');
ylabel('S''(x)');
grid on;

subplot(3, 1, 3);
plot(domain, Szz, 'g', 'LineWidth', 2);
title('Second Derivative: S''''(x)');
xlabel('x');
ylabel('S''''(x)');
grid on;


