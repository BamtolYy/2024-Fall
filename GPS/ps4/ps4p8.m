clear; close all; clc;
%% Get Signal
%----- Setup
Tfull = 0.5;            % Time interval of data to load
fsampIQ = 5.0e6;        % IQ sampling frequency (Hz)
N = floor(fsampIQ*Tfull);
nfft = 2^9;             % Size of FFT used in power spectrum estimation
fIF     = 2.5e6;        % Intermediate frequency (Hz)

%----- Load data
fid = fopen('C:\Users\gsh04\Desktop\2024-Fall\GPS\ps4\niData01head_5MHz.bin','r','l');
Y = fread(fid, [2,N], 'int16')';
Y = Y(:,1) + 1j*Y(:,2);
fclose(fid);
%---- Convert IQ to IF
Tl = 1/fsampIQ;
[Z] = iq2if(real(Y),imag(Y),Tl,fIF);

%% Problem 8

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
Tl = 1/fsampIQ;             % Baseband Sampling time interval in seconds
T  = Tl/2;                  % Bandpass Sampling time interval in seconds
delChip = T/Tc;            % Sampling interval in chips for baseband
Np = 2^nStages - 1;         % Period of the sequence in chips
Ns = length(Z);             % Number of Samples should equal to that of Y(signal)
Ta = 0.001;                 % Accumulation time in seconds
Nk = floor(Ta/T);          % Number of samples in one 1-ms accumulation
% Generate 37 Seqeuences and Oversample them:
codeOS = zeros(Nk,37);
G2tab = [2, 6;3,7;4,8;5,9;1,9;2,10;1,8;2,9;3,10;2,3;3,4;5,6;6,7;7,8;...
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
prn = 2;
% Approximate Doppler (taken from GRID output for PRN 31)
fD = [-1600:10:-1500];
% The Doppler that acquisition and tracking see is opposite fD due to
% high-side mixing
fD_internal = -fD;
% Time vector covering the accumulation
tVec = [0:Nk-1]'*T;
Results = zeros(length(tVec),length(fD_internal));
for m = 1:length(fD_internal)
    for kk = 1:length(tVec)
        jk = round(tVec(kk)*1/T)+1;
        jk
        % Generate the phase argument of the local carrier replica
        ThetaVec = [2*pi*(fIF + fD_internal(m))*tVec];
        % Generate the local carrier replica
        carrierVec = exp(-1i*ThetaVec);
        % Generate the full local replica, with both code and carrier
        lVeck = carrierVec.*codeOS(:,prn);
        % Isolate the kth code interval from the data. xVec here holds the +/-1 and
        % +/-3-valued data samples from dfDataHead.bin.  The first element in xVec
        % holds the first sample in dfDataHead.bin.
        xVeck = Z(jk:jk+Nk-1);
        % Perform correlation and accumulation
        Sk = sum(xVeck.*lVeck);
        % Examine the squared magnitude of Sk in dB.  This should be close to 68.29
        % dB
        SkdB = 10*log10(abs(Sk)^2);
        Results(kk,m) = SkdB;
    end
end
 max(Results(:));

 % for jj=length(Results(1,:))
 % figure,
 % plot(Results(:,jj))
 % end
%--------------------------------------------------------------------------
