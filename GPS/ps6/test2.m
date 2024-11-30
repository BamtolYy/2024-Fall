 close all; clear all; clc;
%% a) Generate 5ms worth of signal and its replica

% ----Genererate Code
% LFSR Parameters:
nStages      = 10;
ciVec1       = [10, 3]';
ciVec2       = [10, 9, 8, 6, 3, 2,]';
a0Vec1       = ones(nStages,1);
a0Vec2       = ones(nStages,1);
% Oversampling Parameters:
fs = 40e6/7;                % Sampling Rate
Tc = 1e-3/1023;             % Chip interval in seconds
T  = 1/fs;                  % Sampling time interval in seconds
delChip = T/Tc;             % Sampling interval in chips
Np = 2^nStages - 1;         % Period of the sequence in chips
Ns = ceil(0.005/T);               % Number of Samples
codeOS = zeros(Ns,1);
G2tab = [2,6];
for j = 1
    [GoldSeq] = generateGoldLfsrSequenceCA(nStages,ciVec1,ciVec2,a0Vec1,...
        a0Vec2,G2tab(j,:));
    % Make code +1/-1 not +1/0
    GoldSeq = 2*GoldSeq - 1;
    % Oversample Code: It makes sense to oversample code, since the code
    % embedded within the signal is sampled at a higher rate than its chip
    % rate. Assuming that the code I generate is sampled at the chip rate,
    % oversampling my code I generated at the rate the signal is sampled
    % will allow my code to correlate with the code embedded in the signal
    GoldSeqOS = oversampleSpreadingCode(GoldSeq,delChip,0,Ns,Np);
    codeOS(:,j) = GoldSeqOS;
end