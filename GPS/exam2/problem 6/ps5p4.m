clear all; close all; clc;
%% White Frequency Noise (Phase Random Walk)
%----Simulation Setup
% Number of simulations
ensemble = 1000;
% Number of Samples per simulation;
Ns = 49999;
% Sampling interval
T = 0.001; % 1ms Given
% Noise parameters
sigma_omega = 0.01;     % radians

%----- Simulation
DeltaOmega = sigma_omega*randn(Ns,ensemble);
DeltaTheta_omega = cumsum(DeltaOmega,1);
for ii= 1: ensemble
    Ccoh(ii) =computeCoherence(DeltaTheta_omega(:,ii),Ns);
end
Ccoh2_mean = mean(Ccoh.^2)
tau = T*Ns;
fprintf(['Coherence Time from white frequency noise: %f \n'],tau)

%% White Frequency Rate Noise (Frequency Random Walk)
%---- Simulation Setup
% Number of simulations
ensemble = 1000;
% Number of Samples per simulation;
Ns_alpha = 1600;
% Sampling interval
T = 0.001; % 1ms Given
% Noise parameters
sigma_alpha = 0.0001;   % radians


%----- Simulation
DeltaAlpha = sigma_alpha*randn(Ns_alpha,ensemble);
DeltaOmega_alpha = cumsum(DeltaAlpha,1);
DeltaTheta_alpha = cumsum(DeltaOmega_alpha,1);
for ii= 1: ensemble
    Ccoh_alpha(ii) =computeCoherence(DeltaTheta_alpha(:,ii),Ns_alpha);
end
Ccoh2_mean_alpha = mean(Ccoh_alpha.^2)
tau_alpha = T*Ns_alpha;
fprintf(['Coherence Time from white frequency Rate noise: %f \n'],tau_alpha)


%% White Phass noise (Phase Random Walk)
%----Simulation Setup
% Number of simulations
ensemble = 1000;
% Number of Samples per simulation;
Ns = 10000;
% Sampling interval
T = 0.001; % 1ms Given
% Noise parameters
sigma = 0.8;     % radians

%----- Simulation
Delta = sigma*randn(Ns,ensemble);
for ii= 1: ensemble
    Ccoh(ii) =computeCoherence(Delta(:,ii),Ns);
end
Ccoh2_mean = mean(Ccoh.^2)
tau = T*Ns;
fprintf(['Coherence Time from white frequency noise: Ccoh^2 stalls at 0.5275'])

fprintf(['\n ---------------------------------------------------\n'])
fprintf(['It doesn''t make sense to estimate the coherence time for the \n'...
        'white phase noise process because with the coherence time, we also seek \n ' ...
        'how many samples we can accumulate before the coherence drops to a value.\n' ...
        'The individual samples of this white noise are completely independent \n' ...
        'of each other and random. The samples are inherently not related to each other.'])

