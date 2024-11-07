close all; clear all; clc;
%% Takeoff Roll

%-------- Parameter Setup
% Atmosphere
rho = 1.18995;                  % kg/m^3
R = 287;                        % J/Kg-K
T =288.15;                      % Kelvin
gamma = 1.4;                     
Papt    = 97716.6;              % Pascal
P0      = 101325;               % Pascal
kdelta  = Papt/P0;              
ft2meter = 0.3048;              % feet to meter conversion rate
h = 532*ft2meter;               % m
ss = sqrt(T*R*gamma);           % m/s
%Airplane
Cd  = 0.025;
Cl  = 0.349;
S   = 125;                      % m^2
% Thrust Paramters
Fstatic = 216000;
K0_h2f = 1;
k1_h2f = 3.281*10^-5;
k2_h2f = 10.764*10^-9;
const_F = kdelta*Fstatic*(K0_h2f+k1_h2f*h+k2_h2f*h^2);
k_0_M2F = 1;
k_1_M2F = -1.07;
k_2_M2F= 0.56;
k_prime_F0 = k_0_M2F;
k_prime_F1 = k_1_M2F/ss;
k_prime_F2 = k_2_M2F/ss^2;
% Other parameters
delt = 0.001;                   % seconds
t= 0:delt:24;
g   = 9.80665;                  % m/s^2
mu  = 0.03;                     % Rolling friction Coefficient
Wdot = 9;                       % N/s

%-------- Variable Setup
p = zeros(length(t),1);
v = zeros(length(t),1);
a = zeros(length(t),1);
W = zeros(length(t),1);
T = zeros(length(t),1);
L = zeros(length(t),1);
N = zeros(length(t),1);
D = zeros(length(t),1);
T(1) = const_F*(k_prime_F0 + k_prime_F1*v(1) + k_prime_F2*v(1)^2);  % N
W(1) = 790100;  % N
N(1) = W(1);
a(1) = g/W(1)*(T(1)-D(1)-mu*N(1));

%-------- Simulate
for i=1:length(t)-1
    p(i+1) = p(i)+v(i)*delt;
    v(i+1) = v(i)+a(i)*delt;
    W(i+1) = W(i)-Wdot*delt;
    T(i+1) = const_F*(k_prime_F0 + k_prime_F1*v(i+1) + k_prime_F2*v(i+1)^2);
    L(i+1) = 0.5*rho*v(i+1)^2*S*Cl;
    N(i+1) = W(i+1)-L(i+1);
    D(i+1) = 0.5*rho*v(i+1)^2*S*Cd;
    a(i+1) = g/W(i+1)*(T(i+1)-D(i+1)-mu*N(i+1));

end

%-------- Plot
subplot(3,1,1)
plot(t,p)
ylabel('Position (m)')
subplot(3,1,2)
plot(t,v)
ylabel('Speed (m/s)')
subplot(3,1,3)
plot(t,a)
ylabel('Acceleration (m/s^2)')
xlabel('Time (s)')
sgtitle('Takeoff Roll')
v(end)
%% Constant Weight Takeoff Roll
%---- Variable Setup
pw = zeros(length(t),1);
vw = zeros(length(t),1);
aw = zeros(length(t),1);
Ww = 790100;                %N
Tw = zeros(length(t),1);
Lw = zeros(length(t),1);
Nw = zeros(length(t),1);
Dw = zeros(length(t),1);

Tw(1) = 216000;  % N

for j = 1:length(t)-1
    pw(j+1) = pw(j)+vw(j)*delt;
    vw(j+1) = g/Ww*(k0*delt-mu*Ww*delt+vw(j)^3/3*(k1+k2-0.5*rho*S*Cd+mu/2*rho*S*Cl));
end

%---- Plot
figure,
subplot(2,1,1)
plot(t,pw)
ylabel('Position (m)')
subplot(2,1,2)
plot(t,vw)
ylabel('Speed (m/s)')
xlabel('Time (s)')
sgtitle('Takeoff Roll Constant Weight')

%% Constant Weight and Acceleration Takeoff Roll
%---- Variable Setup
pa = zeros(length(t),1);
va = zeros(length(t),1);
aa = zeros(length(t),1);
Wa = 790100;                %N

Ta = 216000;  % N
Na = Wa;
a     = g/Wa*(Ta-mu*Na);
for k = 1:length(t)-1
    pa(k+1) = pa(k)+va(k)*delt;
    va(k+1) = va(k)+a*delt;
end

%---- Plot
figure,
subplot(3,1,1)
plot(t,pa)
ylabel('Position (m)')
subplot(3,1,2)
plot(t,va)
ylabel('Speed (m/s)')
subplot(3,1,3)
plot(t,aa)
ylabel('Acceleration (m/s^2)')
xlabel('Time (s)')
sgtitle('Takeoff Roll: Constant Weight and Acceleration')