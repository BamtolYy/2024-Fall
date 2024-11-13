%% Load Data from dfDataHead.bin
% Use the document fftAcqTheory.pdf found on Canvas as your guide.
% Recall that you studied the GP2015 front end in Problem Set 4. The GP2015
% produces digitized data with an intermediate frequency
% fIF = 1.405396825396879 MHz and a sampling rate Ns = 40e6/7 samples per
% second. In the absence of Doppler, there would be Ns/1000 = 40000/7 ≈ 5714 samples per GPS L1 C/A code.
clear; clc;
%----- Setup
Tfull = 0.5;                % Time interval of data to load
fs = 40e6/7;                % Sampling frequency (Hz)
N = fs*Tfull;
N = floor(N/16)*16;         % Number of data samples to load
nfft = 2^10;                % Size of FFT used in power spectrum estimation
fIF  =  1.405396825396879e6; % Hz
%----- Load data
fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps5\dfDataHead.bin"], 'r','l');
[Y,count] = binloadSamples(fid,N,'dual');
Y = Y(:,1);
if(count ~= N)
    error('Insufficient data');
end

%% Genererate Code
%---- Generate all possible PRN (37 SVIDs or PRN Sign No.)
% LFSR Parameters:
nStages      = 10;
ciVec1       = [10, 3]';
ciVec2       = [10, 9, 8, 6, 3, 2,]';
a0Vec1       = ones(nStages,1);
a0Vec2       = ones(nStages,1);
% G2Delay      = [5;6;7;8;17;18;139;140;141;251;252;254;255;256;257;258;...
%     469;470;471;472;473;474;509;512;513;514;515;516;859;860;...
%     861;862;863;950;947;948;950];
% Oversampling Parameters:
Tc = 1e-3/1023;             % Chip interval in seconds
T  = 1/fs;                  % Bandpass Sampling time interval in seconds
delChip = T/Tc;             % Sampling interval in chips
Np = 2^nStages - 1;         % Period of the sequence in chips
Ns = length(Y);             % Number of Samples should equal to that of Y(signal)
Ta = 0.001;                 % Accumulation time in seconds
Nk = floor(Ta/T);           % Number of samples in one 1-ms accumulation
% Generate 37 Seqeuences and Oversample them:
codeOS = zeros(Nk,37);
G2tab = [2,6;3,7;4,8;5,9;1,9;2,10;1,8;2,9;3,10;2,3;3,4;5,6;6,7;7,8;...
    8,9;9,10;1,4;2,5;3,6;4,7;5,8;6,9;1,3;4,6;5,7;6,8;7,9;8,10;1,6;2,7;...
    3,8;4,9;5,10;4,10;1,7;2,8;4,10];
parfor j = 1:length(G2tab)
    [GoldSeq] = generateGoldLfsrSequenceCA(nStages,ciVec1,ciVec2,a0Vec1,...
        a0Vec2,G2tab(j,:));
    % Make code +1/-1 not +1/0
    GoldSeq = 2*GoldSeq - 1;
    % Oversample Code: It makes sense to oversample code, since the code
    % embedded within the signal is sampled at a higher rate than its chip
    % rate. Assuming that the code I generate is sampled at the chip rate,
    % oversampling my code I generated at the rate the signal is sampled
    % will allow my code to correlate with the code embedded in the signal
    GoldSeqOS = oversampleSpreadingCode(GoldSeq,delChip,0,Nk,Np);
    codeOS(:,j) = GoldSeqOS;
end

%--------------------------------------------------------------------------
% % Check if two lfsr m-sequence are truly Golden
% X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
% X2 = generateLfsrSequence(nStages,ciVec2,a0Vec1);
% X1 = 2*X1 - 1;
% X2 = 2*X2 - 1;
% % Calculate t(n)
% if mod(nStages, 2) == 1  % n is odd
%     t_n = 1 + 2^((nStages + 1) / 2);
% else  % n is even
%     t_n = 1 + 2^((nStages + 2) / 2);
% end
% % Define preferred correlation values
% preferredCorrelations = [-t_n, -1, t_n - 2]
% [R12,iiVecSeq] = ccorr(X1,X2);
% figure;clf;
% plot(iiVecSeq,R12);
% title(['Potential X1 and X2 crosscorrelation'])
% ylabel('R_{X1,X2}');
% xlabel('Lag (samples)');
% grid on;
% % Because the crosscorrelation of the two lfsr seqeunce has the expected
% % crosscorrelation values, they do make up gold codes.
%--------------------------------------------------------------------------
prn = 31;
% Approximate Doppler (taken from GRID output for PRN 31)
fD = [-6000:100:6000];
% The Doppler that acquisition and tracking see is opposite fD due to
% high-side mixing
fD_internal = -fD;
% Time vector covering the accumulation
tVec = [0:Nk-1]'*T;
Results = zeros(length(tVec),length(fD_internal));
for m = 1:length(fD_internal)
    for kk = 1:length(tVec)
        jk = round(tVec(kk)*1/T)+1;
        % Generate the phase argument of the local carrier replica
        ThetaVec = [2*pi*(fIF + fD_internal(m))*tVec];
        % Generate the local carrier replica
        carrierVec = exp(-1i*ThetaVec);
        % Generate the full local replica, with both code and carrier
        lVeck = carrierVec.*codeOS(:,prn);
        % Isolate the kth code interval from the data. xVec here holds the +/-1 and
        % +/-3-valued data samples from dfDataHead.bin.  The first element in xVec
        % holds the first sample in dfDataHead.bin.
        xVeck = Y(jk:jk+Nk-1);
        % Perform correlation and accumulation
        Sk = sum(xVeck.*lVeck);
        % Examine the squared magnitude of Sk in dB.  This should be close to 68.29
        % dB
        SkdB = abs(Sk)^2;
        Results(kk,m) = SkdB;
    end
end
figure()
surf(Results)
zlabel('Sk^2')
xlabel('Doppler Frequency, fD, (Hz)')
ylabel('Start Time (s)')
[~,max_index] = max(Results(:));
[ts_index,fD_index]=ind2sub(size(Results),max_index);
apparent_doppler_frequency = fD_internal(fD_index);
start_time = tVec(ts_index+1)*1e6;
sigmaIQ = 130;
CN0 =10*log10((max(Results(:))-2*sigmaIQ^2)/(2*sigmaIQ^2*Ta))
disp(['Apparent Doppler Frequency: ', num2str(apparent_doppler_frequency), ' Hz']);
disp(['Approximate Start Time of First Full C/A Code: ', num2str(start_time), ' microseconds']);
%--------------------------------------------------------------------------
