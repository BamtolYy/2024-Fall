% loadRawSamples
%
% Raw integer samples are read from file byte by byte as
%
% [Stream1_Integer1 Stream2_Integer1 Stream3_Integer1 ...
%  StreamN_Integer1 Stream1_Integer2 Stream2_Integer2 ... ]
%
% Streams are ordered according to the CIRCBUFF_STREAM_IDX configuration
% for each bank; typically the ordering is
% {L1, L1_ALT1, L1_ALT2, ..., L1_ALTM, L2, L2_ALT1, ...}

clear;clc;
%----- Setup
datadir = 'C:\Users\gsh04\Desktop\2024-Fall\GPS\exam3\Problem 4';
filename = 'rawintegersamples_fe.bin';

stream = 1;         % Data stream (between 1 and numStreams)
numStreams = 4;     % Number of data streams
Tfull = 0.5;          % Interval of data to load (sec)
fs = 9.6e6;         % Sampling frequency (Hz)
tSeek = 0;          % Seek time into data (sec)

%----- Load data
fid = fopen([datadir '/' filename], 'r', 'n');
Ns = floor(Tfull*fs);
% numStreams bytes per sample, one for each data stream
seekOffset = floor(tSeek*fs)*numStreams;
status = fseek(fid,seekOffset,-1);
if(status == -1)
  error('tSeek beyond file limit');
end
Y = fread(fid, [numStreams,Ns], 'int8')';
fclose(fid);
if(length(Y(:,1)) < Ns)
  error('Insufficient data');
end
Y = Y(:,stream);


%%  Radiolynx Front-End Characteristics
fIF = 2.391428571429e6 ;   % Hz
T = 1/fs;

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
NC = 10;% Noncoherent sum number
for prn = 8
    % Approximate Doppler (taken from GRID output for PRN 31)
    fD = [0:1/(8*Ta):3000];
    % The Doppler that acquisition and tracking see is opposite fD due to
    % high-side mixing
    fD_internal = fD;
    % Time vector covering the accumulation
    tVec = [0:Nk-1]'*T;
    Sk = zeros(length(tVec),length(fD_internal));
    Sk2 = zeros(length(tVec),length(fD_internal));
    % sigmaIQ2= zeros(length(tVec),length(fD_internal));
    for m = 1:length(fD_internal)
        for kk = 1:length(tVec)
            jk = round(tVec(kk)*1/T)+1;
            % Generate the phase argument of the local carrier replica
            ThetaVec = [2*pi*(fIF + fD_internal(m))*tVec];
            % Generate the local carrier replica
            carrierVec = exp(-1i*ThetaVec);
            % Generate the full local replica, with both code and carrier
            lVeck = carrierVec.*codeOS(:,prn);
            for nn = 1: NC
                % Isolate the kth code interval from the data. xVec here holds the +/-1 and
                % +/-3-valued data samples from dfDataHead.bin.  The first element in xVec
                % holds the first sample in dfDataHead.bin.
                xVeck = Y((nn-1)*Nk+jk:(nn-1)*Nk+jk+Nk-1);
                % Perform correlation and accumulation
                Sk(kk,m) = sum(xVeck.*lVeck);

                % Examine the squared magnitude of Sk in dB.  This should be close to 68.29
                % dB
                Sk2(kk,m) = Sk2(kk,m)+abs(Sk(kk,m))^2;
            end
        end
    end
    figure,
    surf(Sk2)
    zlabel('Sk^2')
    xlabel('Doppler Frequency, fD, (Hz)')
    ylabel('Start Time (s)')
    [~,max_index] = max(Sk2(:));
    [ts_index,fD_index]=ind2sub(size(Sk2),max_index);
    apparent_doppler_frequency = fD(fD_index);
    start_time = tVec(ts_index)*1e6;

    %---- Calculate sigmaIQ^2 from Sk2
    % Define the size of the exclusion region
    region_size = 10;
    % Get the size of the matrix
    [num_rows, num_cols] = size(Sk2);
    % Find the indices of the maximum value
    [row,col] = ind2sub(size(Sk2),max_index);
    % Define the rows and columns to delete
    row_min = max(row - region_size, 1); % Ensure no rows < 1
    row_max = min(row + region_size, num_rows); % Ensure no rows > num_rows
    col_min = max(col - region_size, 1); % Ensure no cols < 1
    col_max = min(col + region_size, num_cols); % Ensure no cols > num_cols
    NoisySk2 = Sk2;
    % Delete the rows and columns
    NoisySk2(row_min:row_max, :) = []; % Remove specified rows
    NoisySk2(:, col_min:col_max) = []; % Remove specified columns
    sigmaIQ2 = mean(NoisySk2(:))/2/NC;

    %---- Calculate C/N0
    CN0 =10*log10((max(Sk2(:))/NC-2*sigmaIQ2)/(2*sigmaIQ2*Ta))


    disp(['Apparent Doppler Frequency: ', num2str(apparent_doppler_frequency), ' Hz']);
    disp(['Approximate Start Time of First Full C/A Code: ', num2str(start_time), ' microseconds']);
    disp(['prn:', num2str(prn)])
    %--------------------------------------------------------------------------
end