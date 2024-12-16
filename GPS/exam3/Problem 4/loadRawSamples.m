% % loadRawSamples
% %
% % Raw integer samples are read from file byte by byte as
% %
% % [Stream1_Integer1 Stream2_Integer1 Stream3_Integer1 ...
% %  StreamN_Integer1 Stream1_Integer2 Stream2_Integer2 ... ]
% %
% % Streams are ordered according to the CIRCBUFF_STREAM_IDX configuration
% % for each bank; typically the ordering is
% % {L1, L1_ALT1, L1_ALT2, ..., L1_ALTM, L2, L2_ALT1, ...}
% 
% clear;clc;
% %----- Setup
% datadir = 'C:\Users\gsh04\Desktop\2024-Fall\GPS\exam3\Problem 4';
% filename = 'rawintegersamples_fe.bin';
% 
% stream = 1;         % Data stream (between 1 and numStreams)
% numStreams = 4;     % Number of data streams
% Tfull = 60;          % Interval of data to load (sec)
% fs = 9.6e6;         % Sampling frequency (Hz)
% tSeek = 0;          % Seek time into data (sec)
% 
% %----- Load data
% fid = fopen([datadir '/' filename], 'r', 'n');
% Ns = floor(Tfull*fs);
% % numStreams bytes per sample, one for each data stream
% seekOffset = floor(tSeek*fs)*numStreams;
% status = fseek(fid,seekOffset,-1);
% if(status == -1)
%   error('tSeek beyond file limit');
% end
% Y = fread(fid, [numStreams,Ns], 'int8')';
% fclose(fid);
% if(length(Y(:,1)) < Ns)
%   error('Insufficient data');
% end
% Y = Y(:,stream);
% 
% 
% %%  Radiolynx Front-End Characteristics
% fIF = 2.391428571429e6 ;   % Hz
% T = 1/fs;
% 
% 
% 
% %% Coarse Search
% disp('----------------------------------------------------------')
% fprintf('                    Coarse Search\n')
% disp('----------------------------------------------------------')
% 
% % Coarse Search Parameter
% fDRange = [-4000:100:4000];
% prn = [1:37];
% NC = 2;                    % Noncoherent sum number
% Ta = 0.01;                 % Accumulation time in seconds
% coarsets = zeros(1, length(prn));
% coarsefD = zeros(1, length(prn));
% % Estimate
% for m = 1:length(prn)
%     [ts, fD,~,~,~] = acquisition(Y,prn(m),fDRange,NC,Ta,fs,fIF);
%     coarsets(m) = ts;
%     coarsefD(m) = fD;
% end

disp('----------------------------------------------------------')

% %% Fine Search
% disp('----------------------------------------------------------')
% fprintf('                    Fine Search\n')
% disp('----------------------------------------------------------')
% 
% prnFine =  find(~isnan(coarsefD));
% Ta =0.01;
% NC = 1;
% tsFine = zeros(length(prnFine),1);
% fDFine = zeros(length(prnFine),1);
% for h = 1:length(prnFine)
%     fDmaxFine = coarsefD(prnFine(h))+100;
%     fDminFine = coarsefD(prnFine(h))-100;
%     fDRangeFine = [fDminFine:1 :fDmaxFine];
%     [ts, fD, Sk2,noiseVariance,sk] = acquisition(Y,prnFine(h),fDRangeFine,NC,Ta,fs,fIF);
%     tsFine(h) = ts;
%     fDFine(h) = fD;
%     peakSk2(h) = Sk2;
%     sigmaIQ2(h) = noiseVariance;
% end

%% (c) Initialize the beat carrier phase estimate
thetaHat = 0;
thetaHat_history(1) = thetaHat;

%% (d) Initialize Moving Window Average

g = 9;                      % PRN index
s.IsqQsqAvg = peakSk2(g);
Sk2avg = s.IsqQsqAvg;
movingAverage = peakSk2(g);
s.Ip = real(sk);
s.Qp = imag(sk);
s.sigmaIQ = sqrt(sigmaIQ2(g));

%% (e)
loopOrder = 3;
Bn_target = 10;
s.Ta = Ta;
s.Tc = 1e-3/1023;             % Chip interval in seconds
[s.Ad,s.Bd,s.Cd,s.Dd,Bn_act] = configureLoopFilter(Bn_target,s.Ta,loopOrder);

%% (f) x_k=0 calculation
vTheta=2*pi*fDFine(g);

[V,D] = eig(s.Ad);
q    = vTheta/(s.Cd*V(:,1));
s.xk = q*V(:,1);

%% (g) Correlate
teml = 0.5;                 % Chips
fc = 1575.42*1e6;
tstart = tsFine(g);
% s.sigmaIQ = sqrt(1.05e5);
Nk = floor(Ta/T);
NumberofAccumulation = floor(length(Y(floor(tsFine(g)/T):end))/Nk);
% vTheta_history = zeros(NumberofAccumulation-1,1);
vTheta_history(1) = vTheta;
Sk2_history    = zeros(NumberofAccumulation-1,1);
Sk2_history(1) = peakSk2(g);

for k = 1 : NumberofAccumulation-1
    [Se_k, Sp_k, Sl_k] = performCorrelations(Y, fs, fIF, tstart, vTheta, thetaHat, teml, prnFine(g), Ta);
    %% (h) Update Moving Window AVerage
    s.Ip = real(Sp_k);
    s.Qp = imag(Sp_k);
    s.Ie = real(Se_k);
    s.Qe = imag(Se_k);
    s.Il = real(Sl_k);
    s.Ql = imag(Sl_k);
    if k <100
        movingAverage(k+1) = abs(Sp_k)^2;
        s.IsqQsqAvg = mean(movingAverage);
        Sk2avg(k+1) = s.IsqQsqAvg/10;
    else
        movingAverage = circshift(movingAverage,-1);
        movingAverage(end) = abs(Sp_k)^2;
        s.IsqQsqAvg = mean(movingAverage);
        Sk2avg(k+1) = s.IsqQsqAvg/10;
    end
    %% (i)
    fc = 1575.42*1e6;
    [s.xk,vTheta] = updatePll(s);
    s.vp = vTheta/(2*pi*fc);
    [vTotal] = updateDll(s);
    vTotal_hist(k+1) = vTotal;
    %% (j) Update Beat Carrier Phase Estimate
    thetaHat = thetaHat+vTheta*Ta;
    vTheta_history(k+1) = vTheta;
    %% (k)
    tstart = tstart -vTotal*Ta+Ta;
    Sk2_history(k+1) = abs(Sp_k)^2;
    Sk_history(k+1) = Sp_k;
    thetaHat_history (k+1) = thetaHat;
    tstart_history(k+1) = tstart;
end
%% i) ii)
figure,
Time = [tsFine(g):Ta:length(vTheta_history)*Ta]';
plot(Time,vTheta_history/(2*pi))
title(['fD for PRN: ', num2str(prnFine(g))]);
% ylim([-2235 -2205])
xlim([0 70])
ylabel ('Hz')
xlabel('Time (sec)')
grid on

%% iii)
figure,
Time = [tsFine(g):Ta:length(vTheta_history)*Ta]';
plot(Time,abs(real(Sk_history)),'black')
hold on,
plot(Time,abs(imag(Sk_history)), 'Color', [.7 .7 .7])
title(['S_k History for PRN: ', num2str(prnFine(g))]);
% ylim([-2235 -2205])
xlim([0 70])
xlabel('Time (sec)')

%% iv) 
CN0 = 10*log10((Sk2avg-2*s.sigmaIQ^2)/(2*s.sigmaIQ^2*Ta));
figure,
plot(Time,CN0)
title(['C/N_0 History for PRN: ', num2str(prnFine(g))]);
ylabel('dB-Hz')
xlabel('Time (sec)')



