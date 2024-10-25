 clear; clc;
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