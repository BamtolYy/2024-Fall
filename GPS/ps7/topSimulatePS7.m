clear; clc;
%% Load Data from dfDataHead.bin
% Use the document fftAcqTheory.pdf found on Canvas as your guide.
% Recall that you studied the GP2015 front end in Problem Set 4. The GP2015
% produces digitized data with an intermediate frequency
% fIF = 1.405396825396879 MHz and a sampling rate Ns = 40e6/7 samples per
% second. In the absence of Doppler, there would be Ns/1000 = 40000/7 ≈ 5714 samples per GPS L1 C/A code.
%----- Setup
Tfull = 10;                % Time interval of data to load
fs = 40e6/7;                % Sampling frequency (Hz)
T = 1/fs;
N = fs*Tfull;
N = floor(N/16)*16;         % Number of data samples to load
nfft = 2^10;                % Size of FFT used in power spectrum estimation
fIF  =  1.405396825396879e6; % Hz
%----- Load data
fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps5\dfDataHead.bin"], 'r','l');
[Y,count] = binloadSamples(fid,N,'dual');

Y = Y(floor(fs*8/16)*16:end,1);
if(count ~= N)
    error('Insufficient data');
end
%% Coarse Search
disp('----------------------------------------------------------')
fprintf('                    Coarse Search\n')
disp('----------------------------------------------------------')

% Coarse Search Parameter
fDRange = [-5000:100:5000];
prn = [1:37];
NC = 10;                    % Noncoherent sum number
Ta = 0.001;                 % Accumulation time in seconds
coarsets = zeros(1, length(prn));
coarsefD = zeros(1, length(prn));
% Estimate
for m = 1:length(prn)
    [ts, fD,~] = acquisition(Y,prn(m),fDRange,NC,Ta,fs,fIF);
    coarsets(m) = ts;
    coarsefD(m) = fD;
end

disp('----------------------------------------------------------')

%% Fine Search
disp('----------------------------------------------------------')
fprintf('                    Fine Search\n')
disp('----------------------------------------------------------')

prnFine =  find(~isnan(coarsefD));
Ta =0.001;
NC = 1;
tsFine = zeros(length(prnFine));
fDFine = zeros(length(prnFine));
for h = 1:length(prnFine)
    fDmaxFine = coarsefD(prnFine(h))+20;
    fDminFine = coarsefD(prnFine(h))-20;
    fDRangeFine = [fDminFine:1:fDmaxFine];
    [ts, fD, Sk2,noiseVariance] = acquisition(Y,prnFine(h),fDRangeFine,NC,Ta,fs,fIF);
    tsFine(h) = ts;
    fDFine(h) = fD;
    peakSk2(h) = Sk2;
    sigmaIQ2(h) = noiseVariance;
end

%% (c) Initialize the beat carrier phase estimate
thetaHat = 0;

%% (d) Initialize Moving Window Average
g = 1;                      % PRN
s.IsqQsqAvg = peakSk2(g);
s.Ip = real(sqrt(peakSk2(g)));
s.Qp = imag(sqrt(peakSk2(g)));
s.sigmaIQ = sqrt(sigmaIQ2(g));
%% (e)
loopOrder = 3;
Bn_target = 10;
s.Ta = Ta;
s.Tc = 1e-3/1023;             % Chip interval in seconds
[s.Ad,s.Bd,s.Cd,s.Dd,Bn_act] = configureLoopFilter(Bn_target,s.Ta,loopOrder);


%% (f) x_k=0 calculation
vk=2*pi*fDFine(g);
[V,D] = eig(s.Ad);
k    = vk/V(:,1)/s.Cd;
s.xk = k*V(:,1);

NumberofAccumulation = Tfull/Ta;
ts = tsFine(g);
for pp = 1 : NumberofAccumulation
    %% (g) Correlate
    teml = 0.5;                 % Chips
    Nk = floor(Ta * fs);        % number of samples in an accumulation
    start = (pp-1)*Nk+1;
    endpt = pp*Nk;
    [Se_k, Sp_k, Sl_k] = performCorrelations(Y(start:endpt), fs, fIF, ts(start), vk, thetaHat(start), teml, g, Ta);

    %% (h) Update Moving Window AVerage
    s.Ip = real(Sp_k);
    s.Qp = imag(Sp_k);
    s.Ie = real(Se_k);
    s.Qe = imag(Se_k);
    s.Il = real(Sl_k);
    s.Ql = imag(Sl_k);

    %% (i)
    [xkp1,s.vp] = updatePll(s);
    [vTotal] = updateDll(s);

    %% (j) Update Beat Carrier Phase Estimate
    if pp == 1
        thetaHat(start:endpt)=thetaHat(start)+vk*[0:T:Ta-T];
    else
        thetaHat(start:endpt)=thetaHat(start-1)+vk*[0:T:Ta-T];
    end
    %% (k)
    ts(start:endpt) = ts(start) -vTotal*[0:T:Ta-T];
end