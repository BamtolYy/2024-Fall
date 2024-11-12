% gridPSD.m
%
% Compute power spectrum of GPS samples taken off the GRID receiver

clear; clc;
%----- Setup
Tfull = 0.5;         % Time interval of data to load
fs = 40e6/7;         % Sampling frequency (Hz)
N = fs*Tfull;        
N = floor(N/16)*16;  % Number of data samples to load
nfft = 2^10;          % Size of FFT used in power spectrum estimation

%----- Load data
fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps5\dfDataHead.bin"], 'r','l'); 
[Y,count] = binloadSamples(fid,N,'dual');
Y = Y(:,1);
if(count ~= N)
  error('Insufficient data');
end