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
N = floor(fs*Tfull);
N = floor(N/16)*16;         % Number of data samples to load
nfft = 2^10;                % Size of FFT used in power spectrum estimation
fIF  =  1.405396825396879e6; % Hz
%----- Load data
% fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps5\dfDataHead.bin"], 'r','l');
fid = fopen(["C:\Users\gsh04\Desktop\2024-Fall\GPS\ps7\dfDataHead.bin"], 'r','l');

[Y,count] = binloadSamples(fid,N,'dual');
acquisitionStartTime = 3;
Y = Y(round(fs*acquisitionStartTime):end,1);
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
Ta =0.01;
NC = 1;
tsFine = zeros(length(prnFine));
fDFine = zeros(length(prnFine));
for h = 1
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
% vTheta=2*pi*-fDFine(g);
% vTheta=2*pi*2200;

vTheta=2*pi*2208;
[V,D] = eig(s.Ad);
q    = vTheta/(s.Cd*V(:,1));
s.xk = q*V(:,1);

%% (g) Correlate
teml = 0.5;                 % Chips
fc = 1575.42*1e6;
tstart = tsFine(g);
% tstart = 5.19e-4;
s.sigmaIQ = sqrt(1.05e5);
Nk = round(Ta/T);
NumberofAccumulation = round(length(Y)/Nk);
vTheta_history = zeros(NumberofAccumulation-1,1);
% vTheta_history(1) = vTheta;
Sk2_history    = zeros(NumberofAccumulation-1,1);
Time = [acquisitionStartTime+tsFine(g):Ta:acquisitionStartTime+length(vTheta_history)*Ta]';
for k = 1 : NumberofAccumulation-1
    if k == 1000
        disp('wait')
    end
    [Se_k, Sp_k, Sl_k] = performCorrelations(Y, fs, fIF, tstart, vTheta, thetaHat, teml, prnFine(g), Ta);
    %% (h) Update Moving Window AVerage
    s.Ip = real(Sp_k);
    s.Qp = imag(Sp_k);
    s.Ie = real(Se_k);
    s.Qe = imag(Se_k);
    s.Il = real(Sl_k);
    s.Ql = imag(Sl_k);
    s.IsqQsqAvg = abs(Sp_k)^2;

    %% (i)
    fc = 1575.42*1e6;
    [xkp1,vTheta] = updatePll(s);
    s.vp = -vTheta/(2*pi*fc);
    [vTotal] = updateDll(s);

    %% (j) Update Beat Carrier Phase Estimate
    thetaHat = thetaHat+vTheta*Ta;
    vTheta_history(k) = vTheta;
    %% (k)
    tstart = tstart -vTotal*Ta+Ta;
    Sk2_history(k) = abs(Sp_k^2);
    Sk_history(k) = Sp_k;
    Skl_history(k) = abs(Sp_k)^2-abs(Sl_k)^2;
    % % Plot IQ plane
    % figure;
    % plot(real(Sk_history(k)), imag(Sk_history(k)), 'o', 'MarkerSize', 5, 'LineWidth', 1.5);
    % hold on;
    % plot([-max(abs(real(Sk_history(k)))) max(abs(real(Sk_history(k))))], [0 0], 'r--'); % Real axis
    % grid on;
    % axis equal;
    % xlabel('Re(S_{p,k})');
    % ylabel('Im(S_{p,k})');
    % title('Phasor Plot of S_{p,k}');

end
% figure,
% plot(Sk2_history)
% title('Sk2')
figure,
plot(-vTheta_history/(2*pi))
title('fD')
ylim([-2235 -2205])
% xlim([0 80])
ylabel ('Hz')
xlabel('Time (sec)')
grid on
% figure,
% plot(Skl_history)
% title('Skp-Skl')

% % Plot IQ plane
% figure;
% plot(real(Sk_history), imag(Sk_history), 'o', 'MarkerSize', 5, 'LineWidth', 1.5);
% hold on;
% plot([-max(abs(real(Sk_history))) max(abs(real(Sk_history)))], [0 0], 'r--'); % Real axis
% grid on;
% axis equal;
% xlabel('Re(S_{p,k})');
% ylabel('Im(S_{p,k})');
% title('Phasor Plot of S_{p,k}');
% hold off;