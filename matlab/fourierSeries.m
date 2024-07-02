function [S, Sz, Szz] = fourierSeries(coefficients, domain, L)
    N = length(coefficients) - 1;
    
    S = zeros(length(domain), 1);     % profile S
    Sz = zeros(length(domain), 1);    % first derivative Sz
    Szz = zeros(length(domain), 1);   % second derivative Szz
    
    for i = 1:length(domain)
        x = domain(i);
        
        % Calculate the series and its derivatives at each point x
        S(i) = coefficients(1);
        
        for n = 1:N
            k = n * pi / L;
            
            S(i) = S(i) + coefficients(n + 1) * cos(k * x);
            Sz(i) = Sz(i) - k * coefficients(n + 1) * sin(k * x);
            Szz(i) = Szz(i) - k^2 * coefficients(n + 1) * cos(k * x);
        end
    end
end

