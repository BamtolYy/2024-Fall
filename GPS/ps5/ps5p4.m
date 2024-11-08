clear all; close all; clc;
%% Simulation Setup
% Number of simulations
ensemble = 4000;
% Number of Samples per simulation;
Ns = 5000;
% Sampling interval
T = 0.001; % 1ms Given
% Noise parameters
sigma_omega = 0.01;     % radians
sigma_alpha = 0.0001;   % radians
%
% Ccoh2_phase    = zeros(Ns,1);
% Ccoh2_freq     = zeros(Ns,1);
% Ccoh2_freqRate = zeros(Ns,1);
DeltaTheta_freq = zeros(Ns, ensemble);
%% Simulation
for m = 1:ensemble
    % Generate white frequency noise
    DeltaOmega = sigma_omega*randn(Ns,1);

    for j = 2:Ns
        DeltaTheta_freq(j,m) = DeltaTheta_freq(j-1,m) + DeltaOmega(j) * T;
    end
    % Generate white frequency rate noise
    % DeltaAlpha = sigma_alpha * randn(Ns, 1);
    % DeltaOmega = zeros(Ns, 1);
    % DeltaTheta_freqRate = zeros(Ns, 1);
    % for j = 2:Ns
    %     DeltaOmega(j) = DeltaOmega(j-1) + DeltaAlpha(j) * T;
    %     DeltaTheta_freqRate(j) = DeltaTheta_freqRate(j-1) + DeltaOmega(j) * T;
    % end
end



% Calculate Ccoh^2
Ccoh2_freq = zeros(Ns,ensemble);
for k = 1:ensemble
    for N = 1:Ns
        Ccoh2_freq(N,k) =  computeCoherence(DeltaTheta_freq(:,k),N)^2;
    end
end


sum1 = sum(Ccoh2_freq,2);
% sum2 =sum(Ccoh2_freqRate,2);

Ccoh2_freq_mean = mean(Ccoh2_freq,2);
% Ccoh2_freqRate_mean = mean(Ccoh2_freqRate,2);
plot(Ccoh2_freq_mean)
xlabel('Number of Samples used for Coherence')
ylabel('Mean C_{coh}^2')
