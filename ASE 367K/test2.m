W = 636.636; %lb
S = 5.500;
b = 195.7;

Ixx = 1.82*10^5;
Izz = 4.97*10^7;
Ixz = 9.70*10^5;

u1 = 210*1.68781;
g = 32.2;
theta_1 = 2.4/180*pi;
alpha_1 = 2.4/180*pi;
rho = 1.2673*10^-3;

Cyb = -0.9;
Cyp = 0;
Cyr = 0;

Cyda = 0;
Cydr = 0.12;

Clb = -0.16;
Clp = -0.34;
Clr = 0.13;

Clda = 0.013;
Cldr = 0.008;

Cnb = 0.16;
Cnp = -0.026;
Cnr = -0.28;

CnTb = 0;
CnTr = 0;

Cnda = 0.0018;
Cndr = -0.1;

Yb = (0.5*rho*u1^2*S/W )*Cyb; 
Yp = (0.5*rho*u1^2*S*b/(2*W*u1))*Cyp;
Yr = (0.5*rho*u1^2*S*b/(2*W*u1))*Cyr;

Yda = (0.5*rho*u1^2*S/W )*Cyda;
Ydr = (0.5*rho*u1^2*S/W )*Cydr;

Lb = (0.5*rho*u1^2*S*b/(Ixx))*Clb;
Lp = (0.5*rho*u1^2*S*b^2/(2*Ixx*u1))*Clp;
Lr =(0.5*rho*u1^2*S*b^2/(2*Ixx*u1))*Clr;

Lda = (0.5*rho*u1^2*S*b/(Ixx))*Clda;
Ldr = (0.5*rho*u1^2*S*b/(Ixx))*Cldr;

Nb = (0.5*rho*u1^2*S*b/(Izz))*Cnb;
Np = (0.5*rho*u1^2*S*b^2/(2*Izz*u1))*Cnp;
Nr = (0.5*rho*u1^2*S*b^2/(2*Izz*u1))*Cnr;

NTb =  (0.5*rho*u1^2*S*b/(Izz))*CnTb;
NTr = (0.5*rho*u1^2*S*b^2/(2*Izz*u1))*CnTr;

Nda = (0.5*rho*u1^2*S*b/(Izz))*Cnda;
Ndr = (0.5*rho*u1^2*S*b/(Izz))*Cndr;


M = [u1 0 0 0 0;0 1 -Ixz/Ixx 0 0; 0 -Ixz/Izz 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
R = [Yb Yp Yr-u1 g*cos(theta_1) 0; Lb Lp Lr 0 0; Nb+NTb Np Nr+NTr 0 0; 0 1 tand(theta_1) 0 0; 0 0 secd(theta_1) 0 0];
F = [Yda Ydr; Lda Ldr; Nda Ndr; 0 0; 0 0];

A = inv(M)*R;
B = inv(M)*F;

C = eye(5);
D = zeros(5,2);