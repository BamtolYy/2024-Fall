% Simulation parameters
numRealizations = 10000; % Number of Monte Carlo realizations
Ns = 1000;               % Number of samples per realization
T = 0.001;               % Sampling interval in seconds (1 ms)
threshold = 0.5;         % Threshold for coherence time estimation

% Noise parameters
sigma_theta = 0.01;      % White phase noise standard deviation (radians)
sigma_omega = 0.01;      % White frequency noise standard deviation (radians per sample interval)
sigma_alpha = 0.0001;    % White frequency rate noise standard deviation (radians per sample interval^2)

% Preallocate arrays for coherence squared values for each noise type
Ccoh2_whitePhase = zeros(Ns, 1);
Ccoh2_whiteFreq = zeros(Ns, 1);
Ccoh2_whiteFreqRate = zeros(Ns, 1);

% Run the simulation for each noise type
for realization = 1:numRealizations
    % 1. White Phase Noise
    DeltaTheta_whitePhase = sigma_theta * randn(Ns, 1);
    for N = 1:Ns
        % Compute coherence for this realization and square it
        Ccoh2_whitePhase(N) = Ccoh2_whitePhase(N) + computeCoherence(DeltaTheta_whitePhase, N)^2;
    end

    % 2. White Frequency Noise (Phase Random Walk)
    DeltaOmega = sigma_omega * randn(Ns, 1);   % Generate white frequency noise
    DeltaTheta_whiteFreq = zeros(Ns, 1);       % Initialize phase
    for j = 2:Ns
        DeltaTheta_whiteFreq(j) = DeltaTheta_whiteFreq(j-1) + DeltaOmega(j) * T; % Euler integration for phase
    end
    for N = 1:Ns
        % Compute coherence for this realization and square it
        Ccoh2_whiteFreq(N) = Ccoh2_whiteFreq(N) + computeCoherence(DeltaTheta_whiteFreq, N)^2;
    end

    % 3. White Frequency Rate Noise (Frequency Random Walk)
    DeltaAlpha = sigma_alpha * randn(Ns, 1);   % Generate white frequency rate noise
    DeltaOmega = zeros(Ns, 1);                 % Initialize frequency
    DeltaTheta_whiteFreqRate = zeros(Ns, 1);   % Initialize phase
    for j = 2:Ns
        DeltaOmega(j) = DeltaOmega(j-1) + DeltaAlpha(j) * T; % Euler integration for frequency
        DeltaTheta_whiteFreqRate(j) = DeltaTheta_whiteFreqRate(j-1) + DeltaOmega(j) * T; % Euler integration for phase
    end
    for N = 1:Ns
        % Compute coherence for this realization and square it
        Ccoh2_whiteFreqRate(N) = Ccoh2_whiteFreqRate(N) + computeCoherence(DeltaTheta_whiteFreqRate, N)^2;
    end
end

% Average coherence squared values over realizations
Ccoh2_whitePhase = Ccoh2_whitePhase / numRealizations;
Ccoh2_whiteFreq = Ccoh2_whiteFreq / numRealizations;
Ccoh2_whiteFreqRate = Ccoh2_whiteFreqRate / numRealizations;

% Estimate coherence time for white frequency noise and white frequency rate noise
Ncoh_freq = find(Ccoh2_whiteFreq <= threshold, 1);
Ncoh_freqRate = find(Ccoh2_whiteFreqRate <= threshold, 1);

Tcoh_freq = Ncoh_freq * T;
Tcoh_freqRate = Ncoh_freqRate * T;

% Display results
fprintf('Estimated coherence time for white frequency noise: %.3f seconds\n', Tcoh_freq);
fprintf('Estimated coherence time for white frequency rate noise: %.3f seconds\n', Tcoh_freqRate);

% Plot coherence function decay
figure;
plot((1:Ns) * T, Ccoh2_whitePhase, 'b', 'DisplayName', 'White Phase Noise');
hold on;
plot((1:Ns) * T, Ccoh2_whiteFreq, 'r', 'DisplayName', 'White Frequency Noise');
plot((1:Ns) * T, Ccoh2_whiteFreqRate, 'g', 'DisplayName', 'White Frequency Rate Noise');
yline(threshold, '--k', 'DisplayName', 'Threshold');
xlabel('Time (seconds)');
ylabel('E[C_{coh}^2(N)]');
legend;
title('Coherence Function Decay for Different Noise Types');
grid on;

