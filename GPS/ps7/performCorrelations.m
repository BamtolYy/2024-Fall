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
% Generate Early, Prompt, and Late C/A Code
[codeEarly, codePrompt, codeLate] = generateEarlyPromptLateCodes(prn, ts, teml, fs, Nk);

% Generate Carrier Replica
ThetaVec = [2*pi*fIF + thetaHat];
carrierVec = exp(-1i*ThetaVec);
lvecEarly  = carrierVec' .* codeEarly;
lvecPrompt = carrierVec' .* codePrompt;
lvecLate   = carrierVec' .* codeLate;
% Perform Correlations
Se_k = sum(Y .* lvecEarly);   % Early
Sp_k = sum(Y .* lvecPrompt); % Prompt
Sl_k = sum(Y .* lvecLate);   % Late
end
