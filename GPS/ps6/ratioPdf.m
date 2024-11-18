function [pw] = ratioPdf(w,sigma1,sigma2,theta1,theta2,ccoef)


a = sqrt(w^2/sigma1^2-2*ccoef*w/(sigma1*sigma2)+1/sigma2^2);
b = theta1*w/sigma1^2-ccoef*(theta1+theta2*w)/(sigma1*sigma2)+theta2/sigma2^2;
c = theta1^2/sigma1^2-2*ccoef*theta1*theta2/(sigma1*sigma2)+theta2^2/sigma2^2;
d = exp((b^2-c*a^2)/(2*(1-ccoef^2)*a^2));
xmin = -b/sqrt((1-ccoef^2)*a);
xmax = b/sqrt((1-ccoef^2)*a);
phiu = @(u) 1/sqrt(2*pi)*exp(-1/2*u.^2);
phiy = integral(phiu,xmin,xmax);
pw = b*d/sqrt(2*pi*sigma1*sigma2*a^3)*phiy+sqrt(1-ccoef^2)/...
    (pi*sigma1*sigma2*a^2)*exp(-c/(2*(1-ccoef^2)));
