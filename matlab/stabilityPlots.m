function [complexPlot, muPlot] = stabilityPlots(lambda)

    % get number of mu 
    Nmu = size(lambda, 2); 

    % Create the domain
    mu = linspace(0.001, 1.0, Nmu);

    % Plot lambda on complex plane
    figure;
    complexPlot = scatter(real(lambda(:)), imag(lambda(:)));
    title('Complex Plane Plot of \lambda');
    xlabel('Real Part');
    ylabel('Imaginary Part');
    legend('off');
    % xlim([0.3, 0.7]);
    % ylim([0.325, 0.35]);

    % Plot max real lambda vs mu
    maxrealLambda = zeros(Nmu, 1);
    parfor i = 1:Nmu
        maxrealLambda(i) = max(real(lambda(:,i)));
    end

    figure;
    muPlot = scatter(mu, maxrealLambda);
    title('Maximum Real Part of \lambda vs \mu');
    xlabel('\mu');
    ylabel('Re{\lambda}');
    legend('off');

    % Combine into one plot (if needed)
    % figure;
    % subplot(1,2,1);
    % scatter(real(lambda(:)), imag(lambda(:)), 1);
    % xlabel('Real Part');
    % ylabel('Imaginary Part');
    % title('Complex Plane Plot of \lambda');
    % subplot(1,2,2);
    % scatter(mu, maxrealLambda, 1);
    % xlabel('\mu');
    % ylabel('Re{\lambda}');
    % title('Maximum Real Part of \lambda vs \mu');
    % sgtitle(['Nmu = ' num2str(Nmu)]);

    return
end
