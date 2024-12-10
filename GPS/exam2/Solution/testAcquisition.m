clear; clc;
%% Load Data from dfDataHead.bin
% Use the document fftAcqTheory.pdf found on Canvas as your guide.
% Recall that you studied the GP2015 front end in Problem Set 4. The GP2015
% produces digitized data with an intermediate frequency
% fIF = 1.405396825396879 MHz and a sampling rate Ns = 40e6/7 samples per
% second. In the absence of Doppler, there would be Ns/1000 = 40000/7 ≈ 5714 samples per GPS L1 C/A code.
%----- Setup
Tfull = 20;                % Time interval of data to load
fs = 40e6/7;                % Sampling frequency (Hz)
T = 1/fs;
N = fs*Tfull;
N = floor(N/16)*16;         % Number of data samples to load
nfft = 2^10;                % Size of FFT used in power spectrum estimation
fIF  =  1.405396825396879e6; % Hz
%----- Load data
fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps5\dfDataHead.bin"], 'r','l');
[Y,count] = binloadSamples(fid,N,'dual');

Y = Y(floor(fs*3/16)*16:end,1);
if(count ~= N)
    error('Insufficient data');
end
% Coherent integration time. CHANGE THIS.
Ta = 1e-2;
% Number of non−coherent integrations. CHANGE THIS.
N = 1;

% ======================================================================= %
%                               IQ to IF
% ======================================================================= %
X = Y;

% ======================================================================= %
%                           Generate Codes
% ======================================================================= %
[prnCodeTable] = generatePrnCodeTable(fs, Ta);

% ======================================================================= %
%                             Acquisition
% ======================================================================= %
% Number of samples in an accumulation.
Nk = floor(Ta * fs);
% Sample times.
tVec = (0 : 1/(fs) : (Nk-1)/(fs))';
% Search vector for Doppler.
dopplerStep = 1/(4*Ta);
dopplerVec = -3.5e3 : dopplerStep :0;
% Search vector for code phase.
codeStep = 3; % Shift by 3 samples
codeVec = 0 : (codeStep / (fs)) : 1e-3 - (1 / (fs)); % All possible start times
for prn = 14
    % Initialize to save time.
    M = zeros(length(codeVec), length(dopplerVec));
    disp(['PRN: ' num2str(prn)]);
    ff = 1; % Frequency bin index
    locCode = prnCodeTable(prn,:)'; % Generate code replica once
    for fD = dopplerVec
        disp(['fD = ' num2str(fD) ' Hz']);
        % Generate local carrier replica
        locExp = exp( -(1i) * (2*pi*fIF*tVec + 2*pi*fD*tVec) );
        % Combined code and carrier local replica
        locReplica = locExp .* locCode;
        % Loop over non−coherent integrations
        for nn = 1 : N
            jj = 1; % code phase index
            for jk = 1 : codeStep : (fs * 1e-3)
                % Dot product for fast performance
                % Add non−coherently
                M(jj, ff) = M(jj, ff) + ...
                    abs( locReplica' * (X(((nn-1)*Nk)+jk : ((nn-1)*Nk)+jk+Nk-1)) )^2;
                jj = jj + 1;
            end
        end
        ff = ff + 1;
    end
    surf(dopplerVec, codeVec, M);
    title(['PRN: ' num2str(prn)]);
    % Ask human to detect acquisition
    acq = input('Is peak detected? (1 = Yes): ');
    if isempty(acq)
        acq = 0;
    end
    if acq
        % Find peak
        [~, maxCol] = max(max(M));
        fDoppler = dopplerVec(maxCol);
        [~, maxRow] = max(max(M'));
        codeStart = codeVec(maxRow);
        % C/N0 computation
        MCrop = M( [1:(maxRow-40), (maxRow+40):end], ...
            [1:(maxCol-(1000/dopplerStep)), (maxCol+(1000/dopplerStep)):end] );
        MCropVec = MCrop(:);
        tsigmaIQ_sq = mean(MCropVec);  % 2*sigma {IQ}ˆ2
        
        C_N0 = 10*log10((max(max(M))-tsigmaIQ_sq)/(tsigmaIQ_sq*Ta));
        
        disp(['Doppler = ' num2str(fDoppler)]);
        disp(['Code Start = ' num2str(codeStart*1e6)]);
        disp(['C/N 0 = ' num2str(C_N0)]);
    end
end