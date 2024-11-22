%% Load Data from rawtrimmed_158.bin
clear; close all; clc;

%% Get Signal
%----- Setup
Tfull = 0.5;         % Time interval of data to load
fs = 40e6/7;         % Sampling frequency (Hz)
N = fs*Tfull;
N = floor(N/16)*16;  % Number of data samples to load
nfft = 2^10;         % Size of FFT used in power spectrum estimation
fIF = 1.610476e6;    % Intermediate frequency (Hz)

%----- Load data
fid = fopen("C:\Users\gsh04\Desktop\2024-Fall\GPS\exam2\problem 5\dataout_raw_trimmed_158.bin", 'r', 'l');
[Y, count] = binloadSamples(fid, N, 'dual');
Y = Y(:, 1);

%% Generate Code
%---- Generate all possible PRN (37 SVIDs or PRN Sign No.)
% LFSR Parameters:
nStages = 10;
ciVec1 = [10, 3]';
ciVec2 = [10, 9, 8, 6, 3, 2]';
a0Vec1 = ones(nStages, 1);
a0Vec2 = ones(nStages, 1);

% Oversampling Parameters:
Tc = 1e-3/1023;          % Chip interval in seconds
T = 1/fs;                % Bandpass Sampling time interval in seconds
delChip = T/Tc;          % Sampling interval in chips
Np = 2^nStages - 1;      % Period of the sequence in chips
Ns = length(Y);          % Number of Samples should equal to that of Y(signal)
Ta = 0.001;              % Accumulation time in seconds
Nk = floor(Ta/T);        % Number of samples in one 1-ms accumulation
NN = 5;                  % Number of code periods for non-coherent accumulation

% Generate 37 Sequences and Oversample them:
codeOS = zeros(Nk, 37);
G2tab = [2, 6; 3, 7; 4, 8; 5, 9; 1, 9; 2, 10; 1, 8; 2, 9; 3, 10; 2, 3; 3, 4; 5, 6; ...
    6, 7; 7, 8; 8, 9; 9, 10; 1, 4; 2, 5; 3, 6; 4, 7; 5, 8; 6, 9; 1, 3; 4, 6; ...
    5, 7; 6, 8; 7, 9; 8, 10; 1, 6; 2, 7; 3, 8; 4, 9; 5, 10; 4, 10; 1, 7; 2, 8; 4, 10];
parfor j = 1:length(G2tab)
    [GoldSeq] = generateGoldLfsrSequenceCA(nStages, ciVec1, ciVec2, a0Vec1, a0Vec2, G2tab(j, :));
    % Make code +1/-1 not +1/0
    GoldSeq = 2 * GoldSeq - 1;
    % Oversample Code
    GoldSeqOS = oversampleSpreadingCode(GoldSeq, delChip, 0, Nk, Np);
    codeOS(:, j) = GoldSeqOS;
end

%% Non-Coherent Accumulation
fD = -300000:1000:0; % Doppler search space
tk = [0:Nk-1]' * T;  % Time vector
sigmaIQ = 149;       % Known noise standard deviation (assumed)
threshold = 39.5;    % Threshold for signal detection
CN0 = zeros(37, 1);  % C/N0 storage for each PRN

for mm = 1:37
    zk2sum = zeros(Nk, 1); % Initialize non-coherent accumulation

    for kk = 1:length(fD)
        Cr = fft(codeOS(:, mm));
        fi = fD(kk) + fIF;

        % Non-coherent accumulation loop
        for ii = 1:NN
            startIdx = 1 + (ii - 1) * Nk;
            endIdx = ii * Nk;

            % Ensure we do not exceed the length of Y
            if endIdx > length(Y)
                break;
            end

            % Perform coherent correlation for each segment
            xkTilde = Y(startIdx:endIdx) .* exp(-1i * 2 * pi * fi * tk);
            XrTilde = fft(xkTilde);
            Zr = XrTilde .* conj(Cr);
            zk = ifft(Zr);
            
            % Accumulate power non-coherently (magnitude squared)
            zk2sum = zk2sum + abs(zk).^2;
        end

        % Find maximum value and its corresponding time index
        [maxValue, kmax] = max(zk2sum);
        CN0(mm) = 10 * log10((maxValue - 2 * sigmaIQ^2) / (2 * sigmaIQ^2 * Ta));

        % Signal detection logic
        if CN0(mm) > threshold
            signalStrength(mm) = CN0(mm);
            start_time(mm) = tk(kmax) * 1e6; % Start time in microseconds
            apparent_fD(mm) = fD(kk);
            
            % Display the results
            disp('----------------------------------------------------------')
            disp(['PRN :', num2str(mm)])
            disp(['Apparent Doppler Frequency: ', num2str(apparent_fD(mm)), ' Hz']);
            disp(['Approximate Start Time from first sample: ', num2str(start_time(mm)), ' microseconds']);
            disp(['C/N0: ', num2str(CN0(mm))])
        end
    end
end
