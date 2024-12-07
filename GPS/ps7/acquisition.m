function [ts,fD] = acquisition(prn,fDRange,NC,Ta)
%% Load Data from dfDataHead.bin
% Use the document fftAcqTheory.pdf found on Canvas as your guide.
% Recall that you studied the GP2015 front end in Problem Set 4. The GP2015
% produces digitized data with an intermediate frequency
% fIF = 1.405396825396879 MHz and a sampling rate Ns = 40e6/7 samples per
% second. In the absence of Doppler, there would be Ns/1000 = 40000/7 ≈ 5714 samples per GPS L1 C/A code.
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

%%
tk = [0:Nk-1]'*T;
threshold = 37;

time = zeros(length(fDRange),1);


for mm = prn
    CN0 = zeros(length(fDRange),1);
    Cr = fft(codeOS(:,mm));
    for kk = 1:length(fDRange)
        zk2sum = zeros(Nk, 1);
        for ii = 1:NC
            jk = (ii-1) * Nk + 1;
            jk_end = ii * Nk;
            fi = -fDRange(kk) + fIF;
            xkTilde = Y(jk:jk_end).*exp(-1i*2*pi*fi*tk);
            XrTilde = fft(xkTilde);
            Zr = XrTilde.*(conj(Cr));
            zk = ifft(Zr);
            zk2sum = zk2sum + abs(zk.^2);
        end
        [maxValue,kmax] = max(zk2sum);
        %---- Calculate sigmaIQ^2 from Sk2
        % Define the size of the exclusion region
        region_size = 100;
        % Define the rows and columns to delete
        row_min = max(kmax - region_size, 1); % Ensure no rows < 1
        row_max = min(kmax + region_size, Nk); % Ensure no rows > num_rows
        NoisyZk2 = zk2sum;
        % Delete the rows and columns
        NoisyZk2(row_min:row_max) = []; % Remove specified rows
        sigmaIQ2 = mean(NoisyZk2(:))/2;
        CN0(kk) = 10*log10((maxValue-2*sigmaIQ2)/(2*sigmaIQ2*Ta));
        time(kk) = kmax;
    end
    [maxCN0,maxfd] = max(CN0(:));
    if  maxCN0 > threshold
        signalStrenghth(mm)=maxCN0;
        start_time(mm) = tk(time(maxfd));
        start_time(mm) = (start_time(mm)-floor(start_time(mm)/(Tc*1023))*Tc*1023)*1e6; % us
        apparent_fD(mm) = fDRange(maxfd);
        disp('----------------------------------------------------------')
        disp(['PRN :',num2str(mm)])
        disp(['Apparent Doppler Frequency: ', num2str(apparent_fD(mm)), ' Hz']);
        disp(['Approximate Start Time from first sample: ', num2str(start_time(mm)), ' microseconds']);
        disp (['C/N0: ', num2str(maxCN0)])
        ts =  start_time(mm);
        fD = apparent_fD(mm);
    else
        ts = NaN;
        fD = NaN;
    end
end

