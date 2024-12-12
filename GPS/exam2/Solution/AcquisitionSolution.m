clear; clc; close all; format long g;
fIF = 5e6;
Tfull = 0.5;
% Time interval of data to load
fsampIQ = 10.0e6;
% IQ sampling frequency (Hz)
N = floor(fsampIQ*Tfull);
nfft = 2^9;
% Size of FFT used in power spectrum estimation
% ======================================================================= %
%                               Load data
% ======================================================================= %
fid = fopen("C:\Users\gsh04\Desktop\2024-Fall\GPS\exam2\Solution\niData03head_10MHz.bin",'r','l');
Y = fread(fid, [2,N], 'int16')';
Y = Y(:,1) + (1i)*Y(:,2);
fclose(fid);
% Coherent integration time. CHANGE THIS.
Ta = 1e-3;
% Number of non−coherent integrations. CHANGE THIS.
N = 1;

% ======================================================================= %
%                               IQ to IF
% ======================================================================= %
X = iq2if(real(Y), imag(Y), 1/fsampIQ, fIF);

% ======================================================================= %
%                           Generate Codes
% ======================================================================= %
[prnCodeTable] = generatePrnCodeTable(2*fsampIQ, Ta);

% ======================================================================= %
%                             Acquisition
% ======================================================================= %
% Number of samples in an accumulation.
Nk = floor(Ta * 2 * fsampIQ);
% Sample times.
tVec = (0 : 1/(2*fsampIQ) : (Nk-1)/(2*fsampIQ))';
% Search vector for Doppler.
dopplerStep = 1/(4*Ta);
dopplerVec = -3.5e3 : dopplerStep : 3.5e3;
% Search vector for code phase.
codeStep = 3; % Shift by 3 samples
codeVec = 0 : (codeStep / (2 * fsampIQ)) : 1e-3 - (1 / (2 * fsampIQ)); % All possible start times
for prn = 3
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
            for jk = 1 : codeStep : (2 * fsampIQ * 1e-3)
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