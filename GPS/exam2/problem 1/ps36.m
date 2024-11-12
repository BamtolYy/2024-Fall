clear all;
close all;
clc;
%% Get Data
load("channel.mat")
M= channel';  % transpose to get in column

%% Ionospheric Delay calc from a model
model = 'broadcast';
ionodata.broadcast.alpha0 = 1.1176e-008;
ionodata.broadcast.alpha1 = 7.4506e-009;
ionodata.broadcast.alpha2 = -5.9605e-008;
ionodata.broadcast.alpha3 = -6.9605e-008;
ionodata.broadcast.beta0  = 90112;
ionodata.broadcast.beta1  = 0;
ionodata.broadcast.beta2  = -296610;
ionodata.broadcast.beta3  = -75536;
tGPS.week    = 1575;
tGPS.seconds = 518201.501;
 rRx = [-742005.851560607;-5462223.38476596; 3198008.7346792]; 

% TXID 
 rSv = [20847329.7083373;-15185642.4780402; 6205281.68907901]; 
[delTauG] = getIonoDelay(ionodata,0,rRx,rSv,tGPS,model);
delTaugG  = round (delTauG,4,'significant');
disp(['--------------------Answer---------------------------------------'])
disp([' delTauG from getIonoDelay: ',num2str(delTauG), ' s'])
% find index of I_L1_p_29 where time equals tGPS.seconds.Note: the seconds 
% do not exactly match each other, so I take only the integer value of 
% GPS.seconds to compare with 4th column of L1_TXID







