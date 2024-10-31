clear all; close all; clc;
m= 636636;
Ixx = 1.82*10^7*32.1740486;
Izz = 4.97*10^7*32.1740486;
Ixz = 9.70*10^5*32.1740486;
S = 5500;
b = 195.7;
u1 = 673.436;
theta = 2.4/180*pi;
rho = 1.2673*10^-3*32.1740486;
q = 1/2*rho*u1^2;
g = 32.17405;

%%
Cyb = -0.9;
Cyp = 0;
Cyr = 0;
Cyda = 0;
Cydr = 0.12;
Clb = -0.16;
Clp = -0.34;
Clr = 0.13;
Clda = 0.013;
Cldr = 0.0080;
Cnb = 0.16;
CnTb = 0;
Cnp = -0.026;
Cnr = -0.28;
CnTr = 0;
Cnda = 0.0018;
Cndr = -0.1;

Yb = q*S/m*Cyb;
Yp = q*S*b/(2*m*u1)*Cyp;
Yr = q*S*b/(2*m*u1)*Cyr;
Yda = q*S/m*Cyda;
Ydr = q*S/m*Cydr;
Lb  = q*S*b/Ixx*Clb;
Lp = q*S*b^2/(2*Ixx*u1)*Clp;
Lr = q*S*b^2/(2*Ixx*u1)*Clr;
Lda = q*S*b/Ixx*Clda;
Ldr = q*S*b/Ixx*Cldr;
Nb  = q*S*b/Izz*Cnb;
NTb = q*S*b/Izz*CnTb;
Np = q*S*b^2/(2*Izz*u1)*Cnp;
Nr = q*S*b^2/(2*Izz*u1)*Cnr;
NTr = q*S*b^2/(2*Izz*u1)*CnTr;
Nda = q*S*b/Izz*Cnda;
Ndr = q*S*b/Izz*Cndr;

M = [u1 0 0 0 0;
    0 1 -Ixz/Ixx 0 0;
    0 -Ixz/Izz 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];

R= [Yb Yp Yr-u1 g*cos(theta) 0;
    Lb Lp Lr 0 0;
    Nb+NTb Np Nr+NTr 0 0;
    0 1 tan(theta) 0 0;
    0 0 sec(theta) 0 0];
F = [Yda Ydr;
    Lda Ldr;
    Nda Ndr;
    0 0;
    0 0];

A = inv(M)*R;
B = inv(M)*F;
% disp(A)
% 
% A = [-0.1067 0 -1 0.0477 0;
%     -2.7427 -0.8404 0.3264 0 0;
%     1.0146 -0.0176 -0.2554 0 0;
%     0 1 0.0419 0 0;
%     0 0 1.0009 0 0];
% B = [0 0.0142;
%     0.2211 0.1482;
%     0.0096 -0.6231;
%     0 0;
%     0 0];
C = eye(5);
D = zeros(5,2);

sys = ss(A,B,C,D);

t = 0:0.01:400;

figure,
impRud(2,:) = zeros(length(t),1);
impRud(1,:) = zeros(length(t),1);
impRud(2,1) = 1/180*pi; % rad
impRudResponse = lsim(sys,impRud,t);
impRudResponse(:,2) = impRudResponse(:,2)*b/(2*u1);
impRudResponse(:,3) = impRudResponse(:,3)*b/(2*u1);
subplot(2,1,1)
plot(t,impRudResponse(:,1:3))
legend('\Delta \beta', '\Delta p/2u1', '\Delta r/2u1s')
ylabel('rad^{-1}')
xlabel('Time (seconds)')
subplot(2,1,2)
plot(t,impRudResponse(:,4:5))
ylabel('rad^{-1}')
xlabel('Time (seconds)')
legend('\Delta \phi', '\Delta \psi')
sgtitle('Rudder Impulse Response')



figure,
impA(2,:) = zeros(length(t),1);
impA(1,:) = zeros(length(t),1);
impA(1,1) = 1/180*pi; % rad
impAResponse = lsim(sys,impA,t);
impAResponse(:,2) = impAResponse(:,2)*b/(2*u1);
impAResponse(:,3) = impAResponse(:,3)*b/(2*u1);
subplot(2,1,1)
plot(t,impAResponse(:,1:3))
legend('\Delta \beta', '\Delta p/2u1', '\Delta r/2u1s')
ylabel('rad^{-1}')
xlabel('Time (seconds)')
subplot(2,1,2)
plot(t,impAResponse(:,4:5))
ylabel('rad^{-1}')
xlabel('Time (seconds)')
legend('\Delta \phi', '\Delta \psi')
sgtitle('Aerlion Impulse Response')


stept = 0:0.01:100;
figure,
stepRud(2,:) = ones(length(stept),1);
stepRud(1,:) = zeros(length(stept),1);
stepRud(2,1:10) = 0;
stepRudResponse = lsim(sys,stepRud,stept);
subplot(2,1,1)
plot(stept,stepRudResponse(:,1:3))
ylabel('rad^{-1}')
xlabel('Time (seconds)')
legend('\Delta \beta', '\Delta p', '\Delta r')
subplot(2,1,2)
plot(stept,stepRudResponse(:,4:5))
ylabel('rad^{-1}')
xlabel('Time (seconds)')
legend('\Delta \phi', '\Delta \psi')
sgtitle('Rudder Step Response')

figure,
stepA(2,:) = zeros(length(stept),1);
stepA(1,:) = ones(length(stept),1);
stepA(1,1:10) = 0;
stepAResponse = lsim(sys,stepA,stept);
subplot(2,1,1)
plot(stept,stepAResponse(:,1:3))
ylabel('rad^{-1}')
xlabel('Time (seconds)')
legend('\Delta \beta', '\Delta p', '\Delta r')
subplot(2,1,2)
plot(stept,stepAResponse(:,4:5))
ylabel('rad^{-1}')
xlabel('Time (seconds)')
legend('\Delta \phi', '\Delta \psi')
sgtitle('Aerlion Step Response')
