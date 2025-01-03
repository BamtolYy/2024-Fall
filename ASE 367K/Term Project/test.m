% function [CP] = findLocCP(V, alpha)
% This functiion determines the center of pressure of TREL rocket, which
% dimensions are provided in the ASE 367K Flight Dynamics Term Project
%
%
% ---- Assumptions:
%       1. Dimensions start from the nose
%       2. At sea level

alpha= 2
V =1
%% Find CP for the individual components
if alpha == 0
    CP = 1;
else
    % CP of Nosecone for Conical Xn
    L = 1.2;                %  m; Length of the nose
    Xn = 0.4666*L;
    Xn = 0.8+6+(1-0.466)*1.2

    % CP of Conical transition XT
    Xp = 7.2;               % m; Nose to d1
    d1 = 0.4;               % m; Forward Diameter
    d2 = 0.25;              % m; Aft Diameter
    LT = 0.8;               % m; d1 to d2
    XT = Xp + LT/3*(1+(1-d1/d2)/(1-(d1/d2)^2));

    % CP of Finset XF
    XB = 6;                 % m; Nose to fin beginning
    cr = 0.8;               % m; Fin root
    ct = 0.7;               % m; Fin tip
    XR = (cr-ct)/2;         % m; Fin root beginning to fin tip beginning
    XF = XR/3*(cr+2*ct)/(cr+ct)+1/6*((cr+ct)-(cr*ct)/(cr+ct));
    XF = 0.8+1.2-XF;

    %% Find normal force coefficient at each component CP
    % Nosecone
    CNaN = 2;

    % Finset
    f= 1;
    R = d1/2;
    dBig = 0.8*sqrt(2);
    s = (dBig - d1)/2;
    theta = 0;
    l = s/cos(theta);
    N = 4;                  % Number of fins
    CNaF = (1+f*R/(s+R))*((4*N*(s/d1)^2)/(1+sqrt(1+(2*l/(cr+ct))^2)));

    % Bonical Boattail
    S1 = pi*d1^2/4;
    S2 = pi*d2^2/4;
    CNaCB = 8/(pi*d1^2)*(S2-S1);


    % %% Find Normal forces at each component CP
    % rho = 1.225; % kg/m^3
    % q = 0.5*rho*V^2;
    % A = pi/4*d1^2;
    % Nc = q*A*alpha*CNaN;
    % Nf = q*A*alpha*CNaF;
    % NCB = q*A*alpha*CNaCB;

    rocket_length = 1.2+6+0.8;
    % CP = (CNaN*Xn + CNaF*XF + NCB*XT)/(Nc+Nf+NCB)
    Cp = (CNaN*(Xn)+CNaF*(XF)+)/(CNaN+CNaF+CNaCB)
end