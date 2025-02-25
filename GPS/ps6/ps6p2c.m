close all; clear all; clc;
%% Parameters
N = 100;
Ta = 0.001;
fs    = 1/Ta;
t = (0:N-1)*Ta';
% True Alpha Parameters
f     = 200;
rho   = 4;
theta = pi/4
nfft = 50*N;



%% Generate Sk
nI = randn(1,N);
nQ = randn(1,N);

k = 0:N-1;
nk = nI +1i*nQ;
Sk = rho*exp(1i*(2*pi*f*k*Ta+theta))+nk;

%% Estimate
% Note the estimator for frequnecy is limited by Ta due to Nyquist
% Frequency

% Frequency: 
Skfft = fft(Sk,nfft)/N;
% freq = fs/nfft*(0:nfft-1);      % Cycles per length of the signal in seconds
freq = fs/nfft*(0:nfft-1); 
% plot(freq,abs(Skfft))
[~,f_ml]=max(abs(Skfft));
disp(['Maximum likelihood estimation of frequency: ', num2str(freq(f_ml))])

% Magnitude:
Skamp = abs(fft(Sk))/N;
% plot(Skamp)
disp(['Maximum likelihood estimation of rho: ', num2str(max(Skamp))])

% Theta:

disp(['Maximum likelihood estimation of theta: ', num2str(angle(Skfft(f_ml)))])
