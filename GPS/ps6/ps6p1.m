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
theta2 = linspace(-pi/2, pi/2, 100); % Avoid exact asymptotes
ftheta = zeros(length(theta2),1);

S = rho*exp(i*theta2);
Q = imag(S);
I = real(S);
W = Q./I;


for ii = 1:length(theta2)
a(ii) = sqrt(W(ii)^2+1);
b(ii) = rho*W(ii)*sin(theta)+rho*cos(theta);
c(ii) = rho^2;
d(ii) = exp((b(ii)^2-c(ii)*a(ii)^2)/(2*a(ii)^2));
phiu = @(u) 1/sqrt(2*pi)*exp(-1/2*u.^2);
xmin = -b(ii)/a(ii);
xmax = b(ii)/a(ii);
phiy = integral(phiu,xmin,xmax);
pw(ii)   = b(ii)*d(ii)/(sqrt(2*pi)*a(ii)^3)*phiy+1/(pi*a(ii)^2)*exp(-c(ii)/2);

ftheta(ii) = pw(ii)*(1+W(ii)^2);
end
figure,
histogram(theta_ML,"Normalization","pdf")
hold on,
plot(theta2,ftheta)