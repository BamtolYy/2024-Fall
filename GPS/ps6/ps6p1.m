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

%% Theoretical
% Assume noise variances are 1. 
% Correlation coefficient between Q and I is zero

w = tan(linspace(-pi/2+0.01,pi/2-0.0,100));
ftheta = zeros(length(w),1);
sigma1 = 1;
sigma2 = 1;
muI = rho*cos(theta);
muQ = rho*sin(theta);
for ii = 1:length(w)
ccoef = 0;
a = sqrt(w(ii)^2/sigma1^2-2*ccoef*w(ii)/(sigma1*sigma2)+1/sigma2^2);
b = muI*(ii)/sigma1^2-ccoef*(muI+muQ*w(ii))/(sigma1*sigma2)+muQ/sigma2^2;
c = muI^2/sigma1^2-2*ccoef*muI*muQ/(sigma1*sigma2)+muQ^2/sigma2^2;
d = exp((b^2-c*a^2)/(2*(1-ccoef^2)*a^2));

muI = rho*cos(theta);
muQ = rho*sin(theta);
xmin = -b/sqrt(1-ccoef^2*a);
xmax = b/sqrt(1-ccoef^2*a);
phi = integral(@(u) exp(-1/2*u.^2)/sqrt(2*pi),xmin,xmax);


% PDF of W=Q/I
pw = b*d/sqrt(2*pi*sigma1*sigma2*a^3)*phi+sqrt(1-ccoef^2)/(pi*sigma1*sigma2*a^2)*exp(-c/(2*(1-ccoef^2)));
ftheta(ii) = pw*abs(1+tan(w(ii))^2);
end
% figure,
% histogram(theta_ML)
% hold on,
plot(w,ftheta)
