clear; close all; clc;

%% Parameters
% Signal Setup
fsampIQ = 5.0e6;          % IQ sampling frequency (Hz)
fIF = 2.5e6;              % Intermediate frequency (Hz)
Ta = 0.001;               % Accumulation interval in seconds (1 ms)
Tl = 1/fsampIQ;           % Sampling interval for baseband samples
T = Tl / 2;               % Bandpass sampling time interval
Nk = floor(Ta / T);       % Number of samples in 1 ms accumulation
Tc = 1e-3 / 1023;         % Chip interval (seconds)

% PRN code parameters
txId = 2;                 % Target PRN number
code_delay_max = Nk;      % Maximum code delay in samples
fD_range = -1600:10:-1500; % Doppler frequency search range (Hz)

%% Load Signal
fid = fopen('C:\Users\gsh04\Desktop\2024-Fall\GPS\ps4\niData01head_5MHz.bin', 'r', 'l');
Y = fread(fid, [2, Nk * 2], 'int16')';
Y = Y(:, 1) + 1j * Y(:, 2);
fclose(fid);
[Z] = iq2if(real(Y), imag(Y), Tl, fIF); % Convert IQ data to IF

%% Generate PRN Code for Target Satellite
nStages = 10;
ciVec1 = [10, 3]';
ciVec2 = [10, 9, 8, 6, 3, 2]';
a0Vec1 = ones(nStages, 1);
a0Vec2 = ones(nStages, 1);
G2tab = [2, 6; 3, 7; 4, 8; 5, 9; 1, 9; 2, 10; 1, 8; 2, 9; 3, 10; 2, 3; ...
         3, 4; 5, 6; 6, 7; 7, 8; 8, 9; 9, 10; 1, 4; 2, 5; 3, 6; 4, 7; 5, 8; ...
         6, 9; 1, 3; 4, 6; 5, 7; 6, 8; 7, 9; 8, 10; 1, 6; 2, 7; 3, 8; 4, 9; ...
         5, 10; 4, 10; 1, 7; 2, 8; 4, 10];
[GoldSeq] = generateGoldLfsrSequenceCA(nStages, ciVec1, ciVec2, a0Vec1, a0Vec2, G2tab(txId, :));
GoldSeq = 2 * GoldSeq - 1; % Convert to ±1
GoldSeqOS = oversampleSpreadingCode(GoldSeq, T/Tc, 0, Nk, 1023); % Oversample the PRN code

%% Initialize Variables for Acquisition
Results = zeros(length(fD_range), code_delay_max);
tVec = (0:Nk-1)' * T; % Time vector for one accumulation interval

%% Acquisition Loop (Doppler and Code Delay Search)
for doppler_idx = 1:length(fD_range)
    fD_internal = -fD_range(doppler_idx); % Account for high-side mixing effect
    ThetaVec = 2 * pi * (fIF + fD_internal) * tVec; % Doppler-compensated phase vector
    carrierVec = exp(-1i * ThetaVec); % Local carrier replica for current Doppler hypothesis
    
    for code_delay = 0:code_delay_max-1
        % Circularly shift the oversampled PRN code by code delay
        local_code = circshift(GoldSeqOS, code_delay);
        
        % Create the full local replica with code and carrier
        local_replica = carrierVec .* local_code;
        
        % Isolate the kth interval from the received signal Z
        xVec = Z(1+code_delay:code_delay+Nk); % Adjust indexing for code delay
        Sk = sum(xVec .* local_replica); % Correlation and accumulation
        
        % Store correlation magnitude
        Results(doppler_idx, code_delay+1) = 10*log10(abs(Sk)^2);
    end
end

%% Identify the Best Doppler and Code Delay Estimates
[max_correlation, max_idx] = max(Results(:));
[best_doppler_idx, best_code_delay_idx] = ind2sub(size(Results), max_idx);

% Estimated Doppler and code delay
estimated_doppler = fD_range(best_doppler_idx);
estimated_code_delay = best_code_delay_idx - 1;

disp(['Estimated Doppler for PRN ', num2str(txId), ': ', num2str(estimated_doppler), ' Hz']);
disp(['Estimated code delay for PRN ', num2str(txId), ': ', num2str(estimated_code_delay), ' samples']);

%% Plot Results
figure;
imagesc(fD_range, 0:code_delay_max-1, Results.');
xlabel('Doppler Shift (Hz)');
ylabel('Code Delay (samples)');
title(['Correlation Magnitude for PRN ', num2str(txId)]);
colorbar;
