syms hx hy hz p q r pd qd rd Ixx Iyy Izz Ixz L M N


IB         = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz];
omegaBdot  = [pd;qd;rd];
omegatilde = [0 -r q; r 0 -p; -q p 0];
omegaB     = [p;q;r];
hB        = [hx; hy; hz];

MB = IB*omegaBdot + omegatilde*(IB*omegaB+hB)

%% HW 5 p3a)
clear all;
close all;
clc; 

cT   = 3.552*10^-5; % kg/(N*s)
S    = 125; % m^2
k    = 0.048;
Cd0  =  0.016;
mi   = 78000; % initial mass kg
mf   = 58000; % final mass kg
rho  = 0.380349568; % air density at 35000 ft. in kg/m^3
g    = 9.7737168; % Gravity Acceleration at 35000 ft. in m/s^2
a    = 295.046; % Speed of sound at 35000 ft. in m/s; constant at 35000 ft
Mach = 0.8; % Given
v    = a*Mach; % V

format longG
fun   = @(m) -1/cT*1./(0.5*rho*v*S*Cd0+2*k*m.^2*g^2./(rho*v^3*S))
range = integral (fun,mi,mf)/1000

% Define the integrand for the range calculation
integrand = @(m) -1 ./ ( (cT / 2) * rho * v * S * ...
                         (Cd0 + k * (2 * (m * g) / (rho * v^2 * S)).^2));
% Calculate the range using the integral
r = integral(integrand, mf, mi)/1000 % From final mass to initial mass
