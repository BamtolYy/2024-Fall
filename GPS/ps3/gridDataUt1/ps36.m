clear all;
close all;
clc;
%% Get Data
load("channel.mat")
M= channel';  % transpose to get in column

% TXID 29------------------------------------------------------------------
iidum  = find(M(:,14)==29 & M(:,10) == 1);
txid29 = M(iidum,:);
% Extract L1 and L2 data
jjdum = find(txid29(:,13) ==0);
kkdum = find(txid29(:,13) == 2);
L1_29 = txid29(jjdum,:);
L2_29 = txid29(kkdum,:);
% Match L1 and L2 that were taken simultaneously and extract their
% pseudoranges and Carrier to Noise ratio
p_29  = []; % Initiate Pseudorange container
c_29  = []; % Initiate Carrier to Noise Ratio container
for i = 1:length(L1_29)
    idx = find( L2_29(:,2) == L1_29(i,2));
    if ~isempty(idx)
        p_29  = [p_29;L1_29(i,4),L1_29(i,8), L2_29(idx,8)]; % Pseudoranges (in meters)
        c_29  = [c_29;L1_29(i,4),L1_29(i,7), L2_29(idx,7)]; % Carrier Phases (in cycles)
    end
end

% TXID 31------------------------------------------------------------------
iidum  = find(M(:,14)==31 & M(:,10) == 1);
txid31 = M(iidum,:);
% Extract L1 and L2 data
lldum = find(txid31(:,13) ==0);
mmdum = find(txid31(:,13) == 2);
L1_31 = txid31(lldum,:);
L2_31 = txid31(mmdum,:);
% Match L1 and L2 that were taken simultaneously and extract their
% pseudoranges and Carrier to Noise ratio
p_31  = []; % Initiate Pseudorange container
c_31  = []; % Initiate Carrier to Noise Ratio container
for i = 1:length(L1_31)
    idx = find( L2_31(:,2) == L1_31(i,2));
    if ~isempty(idx)
        p_31  = [p_31;L1_31(i,4),L1_31(i,8), L2_31(idx,8)]; % Pseudoranges (in meters)
        c_31  = [c_31;L1_31(i,4),L1_31(i,7), L2_31(idx,7)]; % Carrier Phases (in cycles)
    end
end
%% Ionospheric Delay Calc from Measurements

% TXID 29-----------------------------------------------------------------
f_L1 = 1575.42 * 10^6; % Hz
f_L2 = 1227.6 * 10^6; % Hz
% Code Derived Ionospheric Delay
I_L1_p_29 = [p_29(:,1), f_L2^2/(f_L1^2-f_L2^2)*(p_29(:,3)-p_29(:,2))];% Textbook notation of group delay at L1 in meters
% Carrier Derived Ionospheric Delay
% Calculate Wavelengths
lambda_L1      = physconst('LightSpeed') / f_L1; % Wavelength for L1 in meters
lambda_L2      = physconst('LightSpeed') / f_L2; % Wavelength for L2 in meters
I_L1_c_29      = [ c_29(:,1), f_L2^2/(f_L1^2-f_L2^2)*(lambda_L1*c_29(:,2)-lambda_L2*c_29(:,3))];
I_L1_c_29(:,2) = I_L1_c_29(:,2)-I_L1_c_29(1,2); % offset so that it fits with code-derived

% TXID 31------------------------------------------------------------------
% Code Derived Ionospheric Delay
I_L1_p_31 = [p_31(:,1), f_L2^2/(f_L1^2-f_L2^2)*(p_31(:,3)-p_31(:,2))];% Textbook notation of group delay at L1 in meters
% Carrier Derived Ionospheric Delay
% Calculate Wavelengths
lambda_L1      = physconst('LightSpeed') / f_L1; % Wavelength for L1 in meters
lambda_L2      = physconst('LightSpeed') / f_L2; % Wavelength for L2 in meters
I_L1_c_31      = [ c_31(:,1), f_L2^2/(f_L1^2-f_L2^2)*(lambda_L1*c_31(:,2)-lambda_L2*c_31(:,3))];
I_L1_c_31(:,2) = I_L1_c_31(:,2)-I_L1_c_31(1,2); % offset so that it fits with code-derived


%% Plot

% TXID 29-----------------------------------------------------------------
figure(1)
t_29 = L2_29(:,2)-L2_29(1,2);% Create Time vecotr based on measurements from channel.mat
plot(t_29,I_L1_p_29(:,2),t_29,I_L1_c_29(:,2));
legend('Code-Derived','Carrier-Derived ')
hold on,
title('Ionospheric Delay: TXID 29');
xlabel('Time (seconds)');
ylabel('Ionospheric Delay (meters)');

% TXID 31------------------------------------------------------------------
figure(2)
t_31 = L2_31(:,2)-L2_31(1,2);% Create Time vecotr based on measurements from channel.mat
plot(t_31,I_L1_p_31(:,2),t_31,I_L1_c_31(:,2));
legend('Code-Derived','Carrier-Derived ')
hold on,
title('Ionospheric Delay: TXID 31');
xlabel('Time (seconds)');
ylabel('Ionospheric Delay (meters)');
%% Find TEC

% TXID 29-----------------------------------------------------------------
% Pseudorange TEC
TEC_p_29 = [I_L1_p_29(:,1), I_L1_p_29(:,2)*f_L1^2/40.3/10^16]; % in TECU = 10^16 electrons/m^2
% Carrier TEC
TEC_c_29 = [I_L1_c_29(:,1), I_L1_c_29(:,2)*f_L1^2/40.3/10^16]; % in TECU = 10^16 electrons/m^2
figure(3)
plot(t_29,TEC_p_29(:,2))
legend('Code-Derived','Carrier-Derived ')
hold on,
title('Total Electron Content: TXID 29');
xlabel('Time (seconds)');
ylabel('TEC (TECU)');

% TXID 31------------------------------------------------------------------

% Pseudorange TEC
TEC_p_31 = [I_L1_p_31(:,1), I_L1_p_31(:,2)*f_L1^2/40.3/10^16]; % in TECU = 10^16 electrons/m^2
% Carrier TEC
TEC_c_31 = [I_L1_c_31(:,1), I_L1_c_31(:,2)*f_L1^2/40.3/10^16]; % in TECU = 10^16 electrons/m^2
figure(4)
plot(t_31,TEC_p_31(:,2))
legend('Code-Derived','Carrier-Derived ')
hold on,
title('Total Electron Content: TXID 31');
xlabel('Time (seconds)');
ylabel('TEC (TECU)');


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

% TXID 29-----------------------------------------------------------------
rSv_29 = [24597807.6872883; -3065999.1384585; 9611346.77939927]; % in meters ECEF
[delTauG] = getIonoDelay(ionodata,0,rRx,rSv_29,tGPS,model);
ionoDelay_model_29 = delTauG*physconst('LightSpeed') % in meter
% find index of I_L1_p_29 where time equals tGPS.seconds.Note: the seconds 
% do not exactly match each other, so I take only the integer value of 
% GPS.seconds to compare with 4th column of L1_TXID
tdum = find(I_L1_p_29(:,1) == fix(tGPS.seconds)); 
ionoDelay_meas_29   = I_L1_p_29(tdum,2) % meter

% TXID 31-----------------------------------------------------------------
rSv_31 = [2339172.27088689; -16191391.3551878; 21104185.0481546]; % in meters ECEF
[delTauG] = getIonoDelay(ionodata,0,rRx,rSv_31,tGPS,model);
ionoDelay_model_31 = delTauG*physconst('LightSpeed') % in meter
% find index of I_L1_p_31 where time equals tGPS.seconds.Note: the seconds 
% do not exactly match each other, so I take only the integer value of 
% GPS.seconds to compare with 4th column of L1_TXID
tdum = find(I_L1_p_31(:,1) == fix(tGPS.seconds)); 
ionoDelay_meas_31   = I_L1_p_31(tdum,2) % meter




