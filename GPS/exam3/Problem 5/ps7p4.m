clear all; close all; clc;
%%
% Clock Measurements
tR10 = 142874.031415926; % sec
tR18 = 142874.031415926; % sec
tR24 = 142874.031415926; % sec
tR32 = 142874.031415926; % sec

tS10 = 142873.929; % sec
tS18 = 142873.930; % sec
tS24 = 142873.926; % sec
tS32 =  142873.928; % sec

% Pseudorange Calculation Parameters
c = physconst('LightSpeed');
Ip10 = c*1.105596e-08;
Ip18 = c*1.206681e-08;
Ip24 = c*1.582324e-08;
Ip32 = c*1.226210e-08;

T10 = c*1.031992e-08;
T18 = c*1.008367e-08;
T24 = c*1.408656e-08;
T32 = c*1.102516e-08;

dtS10 = 0.000124682;
dtS18 = 0.000324047;
dtS24 = -0.000394024;
dtS32 = -0.000018238;

% Ionospheric and tropospheric delays (meters)
Ip = [Ip10; Ip18; Ip24; Ip32];
T = [T10; T18; T24; T32];

% Satellite clock biases (meters)
dtS = c*[dtS10; dtS18; dtS24; dtS32];

rSvECEF10 =[-7771062.621964, -13392926.713787, 21699866.518681];
rSvECEF18 =[ 4747181.126449, -25802711.537558, 3938318.400306];
rSvECEF24 =[ 14684443.702457, -15351446.405048, 15470725.361304];
rSvECEF32 =[ -15730896.676495, -20849558.582301, 5362760.419279];

rSvECEF = [rSvECEF10; rSvECEF18; rSvECEF24; rSvECEF32];

% % Potential Receiver Locations
% A = [-738641.6553, -5462720.0395, 3197898.1150];
% B = [-741133.4180, -5456322.9775, 3208395.3863];
% C = [-740312.5579, -5457063.2879, 3207249.7207];
% D = [-741991.1305, -5462229.1655, 3198021.3786];
% E = [1101968.2340, -4583484.2540, 4282240.1430];
% F = [6378137.0000, 0.0, 0.0];
% rRx = [A;B;C;D;E;F];
%%
% Caculate Pseudorange Potential Location A

p10 = c*(tR10-tS10);
p18 = c*(tR18-tS18);
p24 = c*(tR24-tS24);
p32 = c*(tR32-tS32);
z =[p10,p18,p24,p32]';



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

x_est(1:3) = round(x_est(1:3),3);
disp('----------------------------------------------------------')

fprintf( 'Pseudoranges')
fprintf(['\n PRN 10: ', num2str(round(p10,3)), ' m'])
fprintf(['\n PRN 18: ', num2str(round(p18,3)), ' m'])
fprintf(['\n PRN 24: ', num2str(round(p24,3)), ' m'])
fprintf(['\n PRN 32: ', num2str(round(p32,3)), ' m','\n'])

disp('----------------------------------------------------------')
fprintf('Reciever ECEF Postion estimate')
fprintf(['\n  x: ',num2str(x_est(1)),'  m'])
fprintf(['\n  y: ',num2str(x_est(2)),' m'])
fprintf(['\n  z: ',num2str(x_est(3)),'  m\n'])
fprintf('%c',948,'t_R : ', round(x_est(4)/c,6), ' seconds')
fprintf('\n')
