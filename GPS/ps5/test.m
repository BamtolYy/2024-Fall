% Parameters
Ns = 1000;               % Number of samples
ensemble = 1000;         % Number of realizations for Monte Carlo
T = 0.001;               % Sampling interval (1 ms)
sigma_theta = 0.01;      % White phase noise standard deviation
sigma_omega = 0.01;      % White frequency noise standard deviation
sigma_alpha = 0.0001;    % White frequency rate noise standard deviation

% Initialize accumulators outside the parfor loop
Ccoh2_phase = zeros(Ns, 1);
Ccoh2_freq = zeros(Ns, 1);
Ccoh2_freqRate = zeros(Ns, 1);

% Simulation using parfor
parfor i = 1:ensemble
    % Temporary variables to accumulate results for each realization
    temp_Ccoh2_phase = zeros(Ns, 1);
    temp_Ccoh2_freq = zeros(Ns, 1);
    temp_Ccoh2_freqRate = zeros(Ns, 1);

    % Generate white phase noise
    DeltaTheta_phase = sigma_theta * randn(Ns, 1);

    % Generate white frequency noise (phase random walk)
    DeltaOmega = sigma_omega * randn(Ns, 1);
    DeltaTheta_freq = zeros(Ns, 1);
    for j = 2:Ns
        DeltaTheta_freq(j) = DeltaTheta_freq(j-1) + DeltaOmega(j) * T;
    end

    % Generate white frequency rate noise (frequency random walk)
    DeltaAlpha = sigma_alpha * randn(Ns, 1);
    DeltaOmega = zeros(Ns, 1);
    DeltaTheta_freqRate = zeros(Ns, 1);
    for j = 2:Ns
        DeltaOmega(j) = DeltaOmega(j-1) + DeltaAlpha(j) * T;
        DeltaTheta_freqRate(j) = DeltaTheta_freqRate(j-1) + DeltaOmega(j) * T;
    end

    % Calculate Ccoh^2 for each N
    for N = 1:Ns
        temp_Ccoh2_phase(N) = computeCoherence(DeltaTheta_phase, N)^2;
        temp_Ccoh2_freq(N) = computeCoherence(DeltaTheta_freq, N)^2;
        temp_Ccoh2_freqRate(N) = computeCoherence(DeltaTheta_freqRate, N)^2;
    end

    % Accumulate the results from this realization
    Ccoh2_phase = Ccoh2_phase + temp_Ccoh2_phase;
    Ccoh2_freq = Ccoh2_freq + temp_Ccoh2_freq;
    Ccoh2_freqRate = Ccoh2_freqRate + temp_Ccoh2_freqRate;
end

% Average over all realizations
Ccoh2_phase = Ccoh2_phase / ensemble;
Ccoh2_freq = Ccoh2_freq / ensemble;
Ccoh2_freqRate = Ccoh2_freqRate / ensemble;
