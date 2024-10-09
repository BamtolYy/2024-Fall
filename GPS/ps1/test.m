% topSimulateTrainDoppler.m
% Top-level script for train Doppler simulation

clear; clc; clf;
%----- Setup
fc = 440;              % Base frequency of the train horn
vTrain = 20;           % Train velocity (m/s)
t0 = 0;                % Initial time
x0 = 0;                % Initial position of the train
delt = 0.01;           % Time step (seconds)
N = 1000;              % Number of samples
vs = 343;              % Speed of sound (m/s)
xObs = 56.5;           % Observer position (m)
dObs = -61;            % Doppler shift (m)

%----- Simulate Doppler Effect
[fDVec, tVec] = simulateTrainDoppler(fc, vTrain, t0, x0, xObs, dObs, delt, N, vs);
fApparentVec = fDVec + fc;  % Apparent frequency

%----- Plot Simulated Data
figure(1)
plot(tVec, fDVec + fc, 'r');
xlabel('Time (seconds)');
ylabel('Apparent horn frequency (Hz)');
grid on;
shg;

%----- Generate Sound Vector
T = delt * N;                     % Total simulation time (sec)
fs = 22050;                       % Sampling frequency (Hz)
deltSamp = 1 / fs;                % Sampling interval (sec)
Ns = floor(T / deltSamp);         % Number of samples
tsamphist = [0:Ns-1]' * deltSamp; % Time history vector
Phihist = zeros(Ns, 1);          % Phase history
fApparentVecInterp = interp1(tVec, fApparentVec, tsamphist, 'spline'); % Interpolate fApparentVec to match time vector

% Generate sound wave based on apparent frequencies
for ii = 2:Ns
    fii = fApparentVecInterp(ii);
    Phihist(ii) = Phihist(ii-1) + 2 * pi * fii * deltSamp;
end
soundVec = sin(Phihist);

%----- Write Sound to Audio File
audiowrite('trainout.wav', soundVec, fs);

%----- Write Frequency-Time History to Output File
save trainData fApparentVec tVec;

%% Load Provided Audio File (H.wav)
%----- Read provided audio file
[y, fs] = audioread('H.wav');    % Read the provided 'H.wav' file

%----- Step 1: Extract Frequency-Time History using Spectrogram -----
% Define spectrogram parameters
windowSize = 2^9;             % Window size for STFT
noverlap = windowSize / 2;      % 50% overlap
nfft = 2^14;  % FFT length (power of 2 for efficiency)

% Compute the spectrogram of the sound data
[S, F, T, P] = spectrogram(y, windowSize, noverlap, nfft, fs);

% Find the frequency with the maximum power at each time point (dominant frequency)
[~, maxIdx] = max(abs(S), [], 1);  % Find index of max power for each time slice
fExtracted = F(maxIdx);            % Extract the corresponding frequency values

%----- Plot the spectrogram and extracted frequency -----
figure;
subplot(2,1,1);
surf(T, F, 10*log10(abs(P)), 'EdgeColor', 'none');
axis tight;
view(0, 90);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of H.wav');

subplot(2,1,2);
plot(T, fExtracted, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Extracted Frequency-Time History from H.wav');

%----- Step 2: Interpolate Extracted Frequency to Match Simulated Data -----
% Interpolate extracted frequency to match the simulation time vector
plotExtracted = interp1(T, fExtracted, tVec, 'spline');

%----- Step 3: Plot Comparison of Simulated and Extracted Data -----
figure;
plot(tVec, plotExtracted, 'b', tVec, fDVec + fc, 'r');
xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
legend('Extracted Frequency from H.wav', 'Simulated Apparent Frequency');
title('Comparison of Extracted and Simulated Frequency-Time Histories');
grid on;
shg;

%----- Step 4: Estimate Observer Position (xObs) and Doppler Shift (dObs) -----
% Cross-correlation to match the extracted and simulated data
fApparentInterp = interp1(tVec, fApparentVec, tVec, 'spline');  % Interpolate simulated data to match extracted time vector

% Cross-correlation between extracted and simulated frequency-time histories
[corr, lag] = xcorr(plotExtracted - mean(plotExtracted), fApparentInterp - mean(fApparentInterp), 'coeff');

% Find the lag that maximizes the cross-correlation (time shift)
[~, idxMaxCorr] = max(corr);
timeShift = lag(idxMaxCorr) * (tVec(2) - tVec(1));  % Time shift corresponding to max correlation

% Estimate the Doppler shift and observer position
% Assuming a rough velocity of the train based on known values (adjust if necessary)
trainVelocity = 20;  % Train velocity in m/s
deltaF = mean(plotExtracted - fApparentInterp);  % Mean frequency shift

% Estimate xObs and dObs based on the Doppler effect
xObsEstimate = (deltaF / (trainVelocity / vs)) * 0.1;  % Precision in 0.1 meters
dObsEstimate = deltaF / fApparentInterp(1);  % Doppler shift estimate

% Display estimated parameters
fprintf('Estimated xObs: %.1f meters\n', xObsEstimate);
fprintf('Estimated dObs: %.1f meters\n', dObsEstimate);

% ----- Play the sound vector
sound([soundVec, y], fs);  % Play both the generated sound and the provided audio

% Plot the provided audio signal
t = linspace(0, length(y) / fs, length(y));  % Time vector for audio file
figure;
plot(t, y);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Original Audio Signal (H.wav)');
grid on;
