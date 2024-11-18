close all; clear all; clc 
%%
ensemble = 1000000;
rho = 4;
theta = 0;
nI = rand(ensemble,1);
nQ = rand(ensemble,1);
I = rho*cos(theta)+nI;
Q = rho*sin(theta)+nQ;
theta_ML = atan2(Q,I);

% Define the parameters
rho = 4; % Signal amplitude
theta = 0; % Initial phase
sigma1 = 1; % Standard deviation for I
sigma2 = 1; % Standard deviation for Q

% Define range of w for the theoretical calculation
w = tan(linspace(-pi/2 + 0.01, pi/2 - 0.01, 10000)); % Avoid singularities at +/- pi/2
ftheta = zeros(length(w), 1);

% Loop through each value of w to calculate the theoretical PDF
for ii = 1:length(w)
    % Calculate a(w), b(w), c, d(w)
    a = sqrt((w(ii)^2 / sigma1^2) + (1 / sigma2^2));
    b = (rho * cos(theta) * w(ii) / sigma1^2) + (rho * sin(theta) / sigma2^2);
    c = (rho^2 * cos(theta)^2 / sigma1^2) + (rho^2 * sin(theta)^2 / sigma2^2);
    d = exp((b^2 - c * a^2) / (2 * a^2));
    
    % Define integration limits for Phi
    xmin = -b / a;
    xmax = b / a;
    
    % Define the integrand for the standard normal PDF
    integrand = @(u) exp(-0.5 * u.^2) / sqrt(2 * pi);
    
    % Numerically integrate to calculate Phi
    phi_pos = integral(integrand, -Inf, xmax); % Integral from -infinity to xmax
    phi_neg = integral(integrand, -Inf, xmin); % Integral from -infinity to xmin
    
    % Calculate Phi difference
    phi = phi_pos - phi_neg;
    
    % Calculate the PDF of W = Q/I using the given formula
    pw = (b * d / (sqrt(2 * pi) * sigma1 * sigma2 * a^3)) * phi + ...
         sqrt(1) / (pi * sigma1 * sigma2 * a^2) * exp(-c / 2);
    
    % Transform to the PDF of theta using the Jacobian of the transformation
    ftheta(ii) = pw * abs(1 + tan(w(ii))^2);
end

%% Plotting the Simulated and Theoretical Distribution
% Simulate the data
ensemble = 1000000;
nI = randn(ensemble, 1); % Noise component for I
nQ = randn(ensemble, 1); % Noise component for Q
I = rho * cos(theta) + nI;
Q = rho * sin(theta) + nQ;
theta_ML = atan2(Q, I);

% Plot the histogram of simulated data
figure;
histogram(theta_ML, 'Normalization', 'pdf', 'BinWidth', 0.01);
hold on;

% Plot the theoretical PDF of theta
plot(linspace(-pi/2 + 0.01, pi/2 - 0.01, 10000), ftheta, 'r', 'LineWidth', 2);

% Set plot labels and legend
xlabel('\theta (radians)');
ylabel('Probability Density');
title('Comparison of Simulated and Theoretical Distributions of \theta_{ML}');
legend('Simulated Histogram', 'Theoretical Distribution');
grid on;
hold off;
