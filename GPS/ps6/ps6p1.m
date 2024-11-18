close all; clear all; clc 
%%
ensemble = 1000000;
rho = 4;
theta = 0;
nI = randn(ensemble,1);
nQ = randn(ensemble,1);
I = rho*cos(theta)+nI;
Q = rho*sin(theta)+nQ;
theta_ML = atan2(Q,I);

%% Theoretical
% % Assume noise variances are 1. 
% % Correlation coefficient between Q and I is zero
% 
sigma1 = 1;
sigma2 = 1;
muI = rho*cos(theta);
muQ = rho*sin(theta);
ccoef = 0;
theta2 = linspace(-pi/2 + 0.01, pi/2 - 0.01, 1000); % Avoid exact asymptotes
ftheta = zeros(length(theta2),1);

for ii = 1:length(theta2)
w(ii)      = tan(theta2(ii));
pw = ratioPdf(w(ii),sigma1,sigma2,muI,muQ,ccoef);
ftheta(ii) = pw*abs(1+tan(theta2(ii))^2);
end


% Plotting the Simulated and Theoretical Distribution
figure,
histogram(theta_ML);
figure,
plot(ftheta)

% Plot the theoretical PDF of theta
% plot(theta2, ftheta, 'r', 'LineWidth', 2);
% 
% % Set plot labels and display the plot
% xlabel('\theta (radians)');
% ylabel('Probability Density');
% title('Comparison of Simulated and Theoretical Distributions of \theta_{ML}');
% legend('Simulated Histogram', 'Theoretical Distribution');
% grid on;
% hold off;