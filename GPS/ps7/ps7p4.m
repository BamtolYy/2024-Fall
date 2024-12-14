clear all; close all; clc;
%%
% Clock Measurements
tR14 = 429615.504253994; % sec
tR21 = 429615.504514099; % sec
tR22 = 429615.505066282; % sec
tR27 = 429615.505190982; % sec

tS14 = 429615.436; % sec
tS21 = 429615.446; % sec
tS22 = 429615.440; % sec
tS27 = 429615.438; % sec

% Pseudorange Calculation Parameters
c = physconst('LightSpeed');
Ip14 = c*2.33838930420718e-08;
Ip21 = c*1.07933295823888e-08;
Ip22 = c*1.69238352966293e-08;
Ip27 = c*1.71386691737267e-08;

T14 = c*1.90304034510942e-08;
T21 = c*8.36092156790743e-09;
T22 = c*1.31971700978366e-08;
T27 = c*1.61167561114859e-08;

dtS14 = 0.000116168688750703;
dtS21 = -0.000109838783947257;
dtS22 = 0.000156846547048633;
dtS27 = 0.000217636932174885;

% Ionospheric and tropospheric delays (meters)
Ip = [Ip14; Ip21; Ip22; Ip27];
T = [T14; T21; T22; T27];

% Satellite clock biases (meters)
dtS = c*[dtS14; dtS21; dtS22; dtS27];

rSvECEF14 =[-15206930.706156, -20832520.6664809, -6469517.7047024];
rSvECEF21 =[ 644129.213511344, -22105331.2073303 14614081.7413721];
rSvECEF22 =[ -18139647.3560073, -12028727.5758169, 15453336.7251709];
rSvECEF27 =[ 17478991.6859777, -16545561.3989088, 11624502.8670108];

rSvECEF = [rSvECEF14; rSvECEF21; rSvECEF22; rSvECEF27];

% Potential Receiver Locations
A = [-738641.6553, -5462720.0395, 3197898.1150];
B = [-741133.4180, -5456322.9775, 3208395.3863];
C = [-740312.5579, -5457063.2879, 3207249.7207];
D = [-741991.1305, -5462229.1655, 3198021.3786];
E = [1101968.2340, -4583484.2540, 4282240.1430];
F = [6378137.0000, 0.0, 0.0];
rRx = [A;B;C;D;E;F];
%%
% Caculate Pseudorange Potential Location A

p14 = c*(tR14-tS14);
p21 = c*(tR21-tS21);
p22 = c*(tR22-tS22);
p27 = c*(tR27-tS27);
z =[p14,p21,p22,p27]';



% Define the nonlinear function h(x)
h = @(x) arrayfun(@(i) ...
    norm(rSvECEF(i, :) - x(1:3)) + x(4) - dtS(i) + Ip(i) + T(i), ...
    1:size(rSvECEF, 1))';
% Initial guess for x = [x, y, z, c*delta_tR]
x0 = [0, 0, 0, 0];

% Solve the system using least squares (nonlinear)
options = optimoptions('lsqnonlin', 'Display', 'iter');
cost_function = @(x) z - h(x); % Residual function
x_est = lsqnonlin(cost_function, x0, [], [], options);

x_est(1:3) = round(x_est(1:3),4);
fprintf(['x: ',num2str(x_est(1)),' y: ', num2str(x_est(2)), ' z: ',num2str(x_est(3))])
sprintf('%c',948,'t_{R} : ', round(x_est(4)/c,6), ' seconds')
