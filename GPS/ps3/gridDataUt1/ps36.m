clear all;
close all;
clc;
%% Get Data
load("channel.mat")
M= channel';  % transpose to get in column
iidum = find(M(:,14)==29 & M(:,10) == 1);
txid29 = M(iidum,:);
% Extract L1 and L2 data
jjdum=find(txid29(:,13) ==0);
kkdum=find(txid29(:,13) == 2);
L1=txid29(jjdum,:);
L2=txid29(kkdum,:);
% Match L1 and L2 that were taken simultaneously and extract their
% pseudoranges and Carrier to Noise ratio
p  = []; % Initiate Pseudorange container
c = []; % Initiate Carrier to Noise Ratio container
for i= 1:length(L1)
    idx = find( L2(:,2) == L1(i,2));
    if ~isempty(idx)
        p  = [p;L1(i,8), L2(idx,8)]; % Pseudoranges (in meters)
        c = [c;L1(i,7), L2(idx,7)]; % Carrier Phases (in cycles)
    end
end

%% Ionospheric Delay calc from a model
model = 'broadcast';
ionodata.broadcast.alpha0 = 4.6566e-009;
ionodata.broadcast.alpha1 = 1.4901e-008;
ionodata.broadcast.alpha2 = -5.9605e-008;
ionodata.broadcast.alpha3 = -5.9605e-008;
ionodata.broadcast.beta0  = 79872;
ionodata.broadcast.beta1  = 65536;
ionodata.broadcast.beta2  = -65536;
ionodata.broadcast.beta3  = -393220;
tGPS.week    = 1490;
tGPS.seconds = 146238.774036515;

rRx    = [1101972.5309609; -4583489.78279095; 4282244.3010423]; % in meters ECEF
rSv_29 = [24597807.6872883; -3065999.1384585; 9611346.77939927]; % in meters ECEF
[delTauG] = getIonoDelay(ionodata,0,rRx,rSv_29,tGPS,model);

%% Ionospheric Delay Calc from Measurements
f_L1 = 1575.42 * 10^6; % Hz
f_L2 = 1227.6 * 10^6; % Hz

% Code Derived Ionospheric Delay
I_L1_p = f_L2^2/(f_L1^2-f_L2^2)*(p(:,2)-p(:,1));% Textbook notation of group delay at L1 in meters

% Carrier Derived Ionospheric Delay
% Calculate Wavelengths
lambda_L1 = physconst('LightSpeed') / f_L1; % Wavelength for L1 in meters
lambda_L2 = physconst('LightSpeed') / f_L2; % Wavelength for L2 in meters
I_L1_c    = f_L2^2/(f_L1^2-f_L2^2)*(lambda_L1*c(:,1)-lambda_L2*c(:,2));
I_L1_c    = I_L1_c-I_L1_c(1); % offset so that it fits with code-derived
%% Plot
figure(1)
t= L2(:,2)-L2(1,2);% Create Time vecotr based on measurements from channel.mat
plot(t,I_L1_p,t,I_L1_c);
legend('Code-Derived','Carrier-Derived ')
hold on,
title('Ionospheric Delay');
xlabel('Time (seconds)');
ylabel('Ionospheric Delay (meters)');

%% Find TEC
% Pseudorange TEC
TEC_p = I_L1_p*f_L1^2/40.3/10^16; % in TECU = 10^16 electrons/m^2
% Carrier TEC
TEC_c = I_L1_c*f_L1^2/40.3/10^16; % in TECU = 10^16 electrons/m^2
figure(2)
plot(t,TEC_p,t,TEC_c)
legend('Code-Derived','Carrier-Derived ')
hold on,
title('Total Electron Content');
xlabel('Time (seconds)');
ylabel('TEC (TECU)');


