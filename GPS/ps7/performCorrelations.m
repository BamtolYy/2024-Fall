function [Se_k, Sp_k, Sl_k] = performCorrelations(Y, fs, fIF, ts, vk, thetaHat, teml, prn, Ta)
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

% Sampling time vector
Nk = floor(Ta * fs); % Number of samples in one accumulation interval
tVec = (0:Nk-1) / fs; % Time vector for accumulation interval

% Generate Early, Prompt, and Late C/A Code
[codeEarly, codePrompt, codeLate] = generateEarlyPromptLateCodes(prn, ts, -teml/2, fs, Nk);

% Generate Carrier Replica
carrierReplica = exp(-1i * (2*pi*fIF*tVec + vk*tVec + thetaHat));

% Perform Correlations
Se_k = sum(Y(1:Nk) .* (carrierReplica .* codeEarly'));   % Early
Sp_k = sum(Y(1:Nk) .* (carrierReplica .* codePrompt')); % Prompt
Sl_k = sum(Y(1:Nk) .* (carrierReplica .* codeLate'));   % Late
end
