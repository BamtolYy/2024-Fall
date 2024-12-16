close all; clear all; clc;
%% Parameters
N = 100;
Ta = 0.001;
fs    = 1/Ta;



load SkSim.mat
N = size(tVec,1);
nfft = 50*N;



% %% Generate Sk
% nI = randn(1,N);
% nQ = randn(1,N);
% 
% k = 0:N-1;
% nk = nI +1i*nQ;
% Sk = rho*exp(1i*(2*pi*f*k*Ta+theta))+nk;

%% Estimate
% Note the estimator for frequnecy is limited by Ta due to Nyquist
% Frequency

% Frequency: 
Skfft = fft(SVec,nfft)/N;
% freq = fs/nfft*(0:nfft-1);      % Cycles per length of the signal in seconds
freq = fs/nfft*(0:nfft-1); 
% plot(freq,abs(Skfft))
[~,f_ml]=max(abs(Skfft));
disp(['Maximum likelihood estimation of frequency: ', num2str(freq(f_ml)), ' Hz'])

% Magnitude:
Skamp = abs(fft(SVec))/N;
% plot(Skamp)
disp(['Maximum likelihood estimation of rho: ', num2str(round(max(Skamp),1))])

% Theta:

disp(['Maximum likelihood estimation of theta: ', num2str(round(angle(Skfft(f_ml)),1)),' radians'])
