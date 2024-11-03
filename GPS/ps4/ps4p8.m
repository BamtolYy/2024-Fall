clear; clc;
%% Get Signal
%----- Setup
Tfull = 0.5;
% Time interval of data to load
fsampIQ = 5.0e6;
% IQ sampling frequency (Hz)
N = floor(fsampIQ*Tfull);
nfft = 2^9;
% Size of FFT used in power spectrum estimation
%----- Load data
fid = fopen('C:\Users\gsh04\Desktop\2024-Fall\GPS\ps4\niData01head_5MHz.bin','r','l');
Y = fread(fid, [2,N], 'int16')';
Y = Y(:,1) + 1j*Y(:,2);
fclose(fid);

%% Problem 8 a) PRN Identifier

% Generate all possible PRN
nStages      = 10;
ciVec1       = [10, 3]';
ciVec2       = [10, 9, 8, 6, 3, 2,]';
a0Vec1       = ones(nStages,1);
a0Vec2       = ones(nStages,1);
G2Delay      = (1:37)';
G2Delay(:,2) = [5;6;7;8;17;18;139;140;141;251;252;254;255;256;257;258;...
    469;470;471;472;473;474;509;512;513;514;515;516;859;860;...
    861;862;863;950;947;948;950];

%--------------------------------------------------------------------------
% Check if two lfsr m-sequence are truly Golden
X1 = generateLfsrSequence(nStages,ciVec1,a0Vec1);
X2 = generateLfsrSequence(nStages,ciVec2,a0Vec1);
X1 = 2*X1 - 1;
X2 = 2*X2 - 1;
% Calculate t(n)
if mod(nStages, 2) == 1  % n is odd
    t_n = 1 + 2^((nStages + 1) / 2);
else  % n is even
    t_n = 1 + 2^((nStages + 2) / 2);
end
% Define preferred correlation values
preferredCorrelations = [-t_n, -1, t_n - 2]
[R12,iiVecSeq] = ccorr(X1,X2);
figure;clf;
plot(iiVecSeq,R12);
title(['Potential X1 and X2 crosscorrelation'])
ylabel('R_{X1,X2}');
xlabel('Lag (samples)');
grid on;
% Because the crosscorrelation of the two lfsr seqeunce has the expected
% crosscorrelation values, they do make up gold codes.
%--------------------------------------------------------------------------

% Generate Gold Sequences for all 37 SVIDs or PRN Sign No.
for i = 1:length(G2Delay)
    [GoldSeq] = generateGoldLfsrSequenceCA(nStages,ciVec1,ciVec2,a0Vec1,...
        a0Vec2,G2Delay(i,2));
    code(:,i) = GoldSeq;
end

for j = 1:length(code)
    [Rseq1,iiVecSeq] = ccorr(Y,X1);
end