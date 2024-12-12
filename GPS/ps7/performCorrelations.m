function [Se_k, Sp_k, Sl_k] = performCorrelations(Y, fs, fIF, ts, vTheta, thetaHat, teml, prn, Ta)
% INPUTS
% Y ---------- Input data samples.
% fs --------- Sampling frequency.
% fIF -------- Intermediate frequency.
% ts --------- Code start time (seconds).
% vk --------- Carrier frequency in radians/sec (v0,k = 2*pi*fD,k).
% thetaHat --- Carrier phase estimate.
% teml ------- Early-minus-late spacing (Chips).
% prn -------- Target PRN.
% Ta --------- Accumulation interval (seconds).
%
% OUTPUTS
% Se_k ------- Early accumulation.
% Sp_k ------- Prompt accumulation.
% Sl_k ------- Late accumulation.
% Number of samples in one accumulation interval
Nk = floor(Ta * fs); 

temlt = teml*1e-3/1023;
jk = floor(ts*fs);
jke = floor((ts-temlt/2)*fs);
jkl = floor((ts+temlt/2)*fs);


% Time vector covering the accumulation
tVecPrompt = [jk+Nk:jk+2*Nk-1]'*1/fs;
tVecEarly = [jke+Nk:jke+2*Nk-1]'*1/fs;
tVecLate = [jkl+Nk:jkl+2*Nk-1]'*1/fs;

% tVecPrompt = [jk:jk+Nk-1]'*1/fs;
% tVecEarly = [jk:jk+Nk-1]'*1/fs;
% tVecLate = [jk:jk+Nk-1]'*1/fs;

% Generate Early, Prompt, and Late C/A Code
[codeEarly, codePrompt, codeLate] = generateEarlyPromptLateCodes(prn, ts, teml, fs, Nk);

% Generate Carrier Replica
ThetaVecPrompt = [2*pi*fIF*tVecPrompt + thetaHat+vTheta*[0:Nk-1]'*1/fs];
ThetaVecEarly = [2*pi*fIF*tVecEarly + thetaHat+vTheta*[0:Nk-1]'*1/fs];
ThetaVecLate = [2*pi*fIF*tVecLate + thetaHat+vTheta*[0:Nk-1]'*1/fs];

% ThetaVecPrompt = [2*pi*fIF*tVecPrompt + vTheta*tVecPrompt];
% ThetaVecEarly = [2*pi*fIF*tVecEarly + vTheta*tVecEarly];
% ThetaVecLate = [2*pi*fIF*tVecLate + vTheta*tVecLate];


carrierVecPrompt = exp(-1i*ThetaVecPrompt);
carrierVecEarly = exp(-1i*ThetaVecEarly);
carrierVecLate = exp(-1i*ThetaVecLate);

lvecPrompt = carrierVecPrompt .* codePrompt;
lvecEarly  = carrierVecEarly .* codeEarly;
lvecLate   = carrierVecLate .* codeLate;

% Perform Correlations
xVecPrompt = Y(jk:jk+Nk-1);
xVecEarly = Y(jke:jke+Nk-1);
xVecLate = Y(jkl:jkl+Nk-1);

Sp_k = sum(xVecPrompt .* lvecPrompt); % Prompt
Se_k = sum(xVecEarly .* lvecEarly);   % Early
Sl_k = sum(xVecLate .* lvecLate);   % Late
% disp(['prompt'])
abs(Sp_k)^2;
abs(Se_k)^2;
abs(Sl_k)^2;

end
