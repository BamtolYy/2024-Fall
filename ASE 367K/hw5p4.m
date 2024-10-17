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

cT   = 3.522*10^-5; % kg/(N*s)
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

fun   = @(m) 1./(0.5*rho*v*S*Cd0+(k*m.^2*g^2)/(2*rho*v^3*S))*1./-cT;
range = integral (fun,mi,mf)


