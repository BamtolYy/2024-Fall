clear all; close all; clc;
u1    = 673.436; % ft/s
g     = 32.2;  % ft/s^2
S     = 5500;  % ft^2
b     = 195.7; % ft
m     = 636636; % lb
theta = 2.4/180*pi; % rad
Ixz = 9.70*10^5*g; % lbs*ft^2
Ixx = 1.82*10^7*g; % lbs*ft^2
Izz = 4.97*10^7*g; % lbs*ft^2
rho = 0.0407741717; % lb/ft^3
q1 = rho*u1^2/2;


%% Dimensional stability Derivatives
Cyb = -0.9;
Cyp = 0; 
Cyr = 0;
Cyda = 0;
Cydr = 0.12;
Clb  = -0.16;
Clp  = -0.34;
Clr = 0.13;
Clda  = 0.013;
Cldr = 0.008;
Cnb = 0.16;
CnTb = 0;
Cnp = -0.026;
Cnr = -0.28;
CnTr = 0;
Cnda = 0.0018;
Cndr = -0.10;

Yb  = q1*S/m*Cyb;
Yp  = q1*S*b/(2*m*u1)*Cyp;
Yr  = q1*S*b/(2*m*u1)*Cyr;
Yda = q1*S/m*Cyda;
Ydr = q1*S/m*Cydr;
Lb  = q1*S*b/Ixx*Clb;
Lp  = q1*S*b^2/(2*Ixx*u1)*Clp;
Lr  = q1*S*b^2/(2*Ixx*u1)*Clr;
Lda = q1*S*b/Ixx*Clda;
Ldr = q1*S*b/Ixx*Cldr;
Nb  = q1*S*b/Izz*Cnb;
NTb = q1*S*b/Izz*CnTb;
Np  = q1*S*b^2/(2*Izz*u1)*Cnp;
Nr  = q1*S*b^2/(2*Izz*u1)*Cnr;
NTr = q1*S*b^2/(2*Izz*u1)*CnTr;
Nda = q1*S*b/Izz*Cnda;
Ndr = q1*S*b/Izz*Cndr;

%% Set up state space equation

M = [u1 0 0 0 0;
    0 1 -(Ixz/Ixx) 0 0;
    0 -(Ixz/Izz) 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];
R = [Yb Yp Yr-u1 g*cos(theta) 0;
    Lb Lp Lr 0 0;
    Nb+NTb Np Nr+NTr 0 0;
    0 1 tan(theta) 0 0;
    0 0 sec(theta) 0 0];
F= [Yda Ydr;
    Lda Ldr;
    Nda Ndr;
    0 0;
    0 0];

A = inv(M)*R;
B = inv(M)*F;
C = eye(5);
D = zeros(5,1);

%% Problem 1 a)
[eVec,eVal] = eig(A);
% Make them unitless
eVecUnitless = eVec;
eVecUnitless(2,:) = abs(eVecUnitless(2,:)*b/(2*u1)); % Not sure how we get this factor?
eVecUnitless(3,:) = abs(eVecUnitless(3,:)*b/(2*u1)); % Not sure how we get this factor?

eVecUnitless(:,2) = eVecUnitless(:,2)/eVecUnitless(4,2);
eVecUnitless(:,3) = eVecUnitless(:,3)/eVecUnitless(5,3);
eVecUnitless(:,4) = eVecUnitless(:,4)/eVecUnitless(4,4);
eVecUnitless(:,5) = eVecUnitless(:,5)/eVecUnitless(4,4);

% iidum = eVecUnitless(2,2)/0.1364;
% jjdum = eVecUnitless(3,3)/0.0032;
% eVecUnitless(2,:) = abs(eVecUnitless(2,:)/iidum); % Not sure how we get this factor?
% eVecUnitless(3,:) = abs(eVecUnitless(3,:)/jjdum); % Not sure how we get this factor?

disp(eVecUnitless)
% disp(iidum)
% disp(jjdum)
% Normalize respect to delta theta
% 
% eVec12 = abs(eVecUnitless(:,1)/eVecUnitless(4,1));
% eVec34 = abs(eVecUnitless(:,3)/eVecUnitless(4,3));
