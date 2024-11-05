clear all; close all; clc;
%% Simulation Setup
% Number of simulations
ensemble = 4000;
% Number of Samples per simulation;
Ns = 10000;
% Sampling interval
T = 0.001; % 1ms Given
% Noise parameters
sigma_theta = 0.01;
sigma_omega = 0.01;
sigma_alpha = 0.0001;
%
Ccoh2_phase    = zeros(Ns,1);
Ccoh2_freq     = zeros(Ns,1);
Ccoh2_freqRate = zeros(Ns,1);

%% Simulation
for i = 1:ensemble
    % Generate white phase noise
    DeltaTheta_phase = sigma_theta*rand(Ns,1);
    % Generate white frequency noise
    DeltaOmega = sigma_omega*rand(Ns,1);
    DeltaTheta_freq = zeros(Ns, 1);
    for j = 2:Ns
        DeltaTheta_freq(j) = DeltaTheta_freq(j-1) + DeltaOmega(j) * T;
    end
    % Generate white frequency rate noise
    DeltaAlpha = sigma_alpha * randn(Ns, 1);
    DeltaOmega = zeros(Ns, 1);
    DeltaTheta_freqRate = zeros(Ns, 1);
    for j = 2:Ns
        DeltaOmega(j) = DeltaOmega(j-1) + DeltaAlpha(j) * T;
        DeltaTheta_freqRate(j) = DeltaTheta_freqRate(j-1) + DeltaOmega(j) * T;
    end

    % Calculate Ccoh^2
    for N = 1:Ns
        Ccoh2_phase(N) = Ccoh2_phase(N)+ computeCoherence(DeltaTheta_phase,N)^2;
        Ccoh2_freq(N) = Ccoh2_freq(N)+ computeCoherence(DeltaTheta_freq,N)^2;
        Ccoh2_freqRate(N) = Ccoh2_freqRate(N)+ computeCoherence(DeltaTheta_freqRate,N)^2;
    end

end

