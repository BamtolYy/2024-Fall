clear; clc;
%%
%----- Setup
Tfull = 0.5;           % Time interval of data to load
fsampIQ = 5.0e6;       % IQ sampling frequency (Hz)


N = floor(fsampIQ*Tfull);
nfft = 2^9;            % Size of FFT used in power spectrum estimation
%----- Load data
fid = fopen('C:\Users\gsh04\Desktop\2024-Fall\GPS\ps4\niData01head_5MHz.bin','r','l');
Y = fread(fid, [2,N], 'int16')';
Y = Y(:,1) + j*Y(:,2);
fclose(fid);
fs = 4/0.8;
%%
%----- Compute power spectrum estimate
[Syy,fVec] = pwelch(Y,hann(nfft),[],nfft,fs);

%%
%----- Plot results
yLow = -180;
yHigh = -130;
T = nfft/fs;
delf = 1/T;
fcenter = (nfft/2)*delf;
fVec = fVec - fcenter;
Syy = [Syy(nfft/2 + 1 : end); Syy(1:nfft/2)];
area(fVec/1e6,10*log10(Syy),yLow);
% ylim([yLow,yHigh]);
grid on;
shg;
xlabel('Frequency (MHz)');
ylabel('Power density (dB/Hz)');
figset
title('Power spectral density estimate of GPS L1 Signal');
shg;
