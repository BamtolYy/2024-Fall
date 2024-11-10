clear all; close all; clc;

%% Simulation Setup
% Number of samples
Ns = 1000;
T = 0.001; % 1 ms sampling interval
sigma_omega = 0.01; % radians per sampling interval

% Small drift term to simulate natural phase drift over time
phase_drift = 0.0001 * (1:Ns)';

% Generate white frequency noise and add drift
DeltaTheta_freq = phase_drift + sigma_omega * randn(Ns, 1) * T;

% Preallocate coherence result
Ccoh2_freq = zeros(Ns, 1);

% Calculate coherence for each N up to Ns
for N = 1:Ns
    % Compute coherence
    Ccoh2_freq(N) = abs((1/N) * sum(exp(1i * DeltaTheta_freq(1:N))))^2;
end

% Plot the coherence function
figure;
plot(Ccoh2_freq)
xlabel('Number of Samples used for Coherence')
ylabel('C_{coh}^2')
title('Coherence with Small Phase Drift')
grid on;

% Plot complex exponentials to visualize phase variation
figure;
plot(real(exp(1i * DeltaTheta_freq)), imag(exp(1i * DeltaTheta_freq)), '.')
xlabel('Real Part')
ylabel('Imaginary Part')
title('Complex Exponentials of Phase Terms')
axis equal;
grid on;
