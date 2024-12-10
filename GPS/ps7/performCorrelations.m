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

temlt = teml*1e-3/1023;
jk = round(ts*fs);
jke = round((ts-temlt)*fs);
jkl = round((ts+temlt)*fs);
% Number of samples in one accumulation interval
Nk = round(Ta * fs); 

% Time vector covering the accumulation
tVec = [jk:jk+Nk-1]'*1/fs;


% Generate Early, Prompt, and Late C/A Code
[codeEarly, codePrompt, codeLate] = generateEarlyPromptLateCodes(prn, ts, teml, fs, Nk);

% Generate Carrier Replica
ThetaVec = [2*pi*fIF*tVec + thetaHat+vTheta*[0:Nk-1]'*1/fs];
carrierVec = exp(-1i*ThetaVec);
lvecEarly  = carrierVec .* codeEarly;
lvecPrompt = carrierVec .* codePrompt;
lvecLate   = carrierVec .* codeLate;
% Perform Correlations
xVecPrompt = Y(jk:jk+Nk-1);
xVecEarly = Y(jke:jke+Nk-1);
xVecLate = Y(jkl:jkl+Nk-1);
Se_k = sum(xVecEarly .* lvecEarly);   % Early
Sp_k = sum(xVecPrompt .* lvecPrompt); % Prompt
Sl_k = sum(xVecLate .* lvecLate);   % Late
% disp(['prompt'])
abs(Sp_k)^2
abs(Se_k)^2
abs(Sl_k)^2

end
