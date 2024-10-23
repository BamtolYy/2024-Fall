Xu     =-0.0188;
XTu    = 0;
Xa     = 11.5905;
Xde    = 0;
Zu     = -0.1862;
Za     = -149.4408;
Zq     = -6.8045;
Zadot  = -8.4426;
Zde    = -8.7058;
Mu     = 0.0001;
MTu    = 0;
Ma     = -0.5294;
MTa    = 0;
Mq     = -0.4275;
Madot  = -0.0658;
Mde    = -0.5630;
u1     = 279.1094; % ft/s
theta1 = 0;
g = 32.2; % ft/s^2
%% Set Up Matrix A, B C, D

R = [Xu+XTu Xa 0 -g*cos(theta1);...
    Zu Za u1+Zq -g*sin(theta1);...
    Mu+MTu Ma+MTa Mq 0;...
    0 0 1 0];

M = [1 0 0 0;
    0 u1-Zadot 0 0;...
    0 -Madot 1 0;...
    0 0 0 1];

F = [Xde; Zde; Mde; 0];

A = inv(M)*R;
B = inv(M)*F;

C = eye(4);

D = zeros(4,1);


%% eig Vectors and Values
[eVec,eVal] = eig(A);

% Make them unitless
eVecUnitless = eVec;
eVecUnitless(1,:) = eVec(1,:)/u1;
eVecUnitless(3,:) = eVec(3,:)/20; % Not sure how we get 20?
% eVecUnitless

% Normalize respect to delta theta

eVec12 = abs(eVecUnitless(:,1)/eVecUnitless(4,1)); 
eVec34 = abs(eVecUnitless(:,3)/eVecUnitless(4,3)); 


%% Answer Questions

% Short
wd   = imag(eVal(1,1));
wn   = sqrt((real(eVal(1,1)))^2+(wd)^2);
damp = abs(real(eVal(1,1))/wn);
delT = abs(log(2)/(real(eVal(1,1))));
N = abs(log(2)*wd/(2*pi*real(eVal(1,1))));
% charEqShort=conv([1 -eVal(1,1)],[1 -eVal(2,2)]);
% wn   = sqrt(charEqShort(end));
% damp = charEqShort(2)/2/wn;
% wd   = wn*sqrt(1-damp^2);
% beta = atan(sqrt(1-damp^2)/damp); 
% tr   = (pi-beta)/wd; % Rise time (0-100%)
% tp   = pi/wd;% Peak Time
% Mp   = exp(-damp*pi/sqrt(1-damp^2)); % Maximum overshoot
% ts   = 4/(damp*wn); % Settling Time (2% criterion)
% delT = log(2)/(damp*wn);
% N = log(2)*wd/(2*pi*damp*wn);

% Long
wdLong   = imag(eVal(3,3));
wnLong   = sqrt((real(eVal(3,3)))^2+(wdLong)^2);
dampLong = abs(real(eVal(3,3))/wnLong);
delTLong = abs(log(2)/(real(eVal(3,3))));
NLong = abs(log(2)*wdLong/(2*pi*real(eVal(3,3))));
%% b)
[b,a] = ss2tf(A,B,C,D);


